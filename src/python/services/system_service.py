"""
System resource monitoring service.

This service encapsulates system resource monitoring and diagnostic logic,
providing a unified interface for both API endpoints and AI agent tools.
"""

import logging
import os
import sys
import json
import multiprocessing
import re
import importlib
import shutil
import site
import subprocess
from datetime import datetime
from typing import Dict, Any, Optional, Tuple, List
from importlib import metadata, util

from quantum_calc import get_process_manager, get_current_settings
from quantum_calc.resource_manager import get_resource_manager
from quantum_calc._calculation_repository import CalculationRepository
from .exceptions import ServiceError, ValidationError

logger = logging.getLogger(__name__)

PIP_INSTALL_TIMEOUT_SECONDS = 30 * 60
CUDA_DETECTION_TIMEOUT_SECONDS = 5

CUDA_PACKAGE_MAP = {
    11: ("gpu4pyscf-cuda11x", "cutensor-cu11"),
    12: ("gpu4pyscf-cuda12x", "cutensor-cu12"),
    13: ("gpu4pyscf-cuda13x", "cutensor-cu13"),
}
CUPY_CUTENSOR_RECOMMENDED_BY_CUDA = {
    11: [
        ("13.4.1", "2.2.0"),
        ("13.3.0", "2.0.2"),
    ],
    12: [
        ("13.4.1", "2.2.0"),
        ("13.3.0", "2.0.2"),
    ],
}
LIBXC_VERSION_BY_CUDA = {
    11: "0.5",
    12: "0.5",
    13: "0.7",
}


class SystemService:
    """Service for system resource monitoring and diagnostics."""

    def _detect_cuda_version(self) -> Tuple[Optional[str], Optional[int], Optional[int], Optional[str]]:
        if shutil.which("nvcc") is None:
            return None, None, None, "nvcc command not found"

        try:
            result = subprocess.run(
                ["nvcc", "--version"],
                capture_output=True,
                text=True,
                timeout=CUDA_DETECTION_TIMEOUT_SECONDS
            )
            output = (result.stdout or "") + "\n" + (result.stderr or "")
        except subprocess.TimeoutExpired:
            return (
                None,
                None,
                None,
                f"nvcc --version timed out after {CUDA_DETECTION_TIMEOUT_SECONDS} seconds"
            )

        if result.returncode != 0:
            message = (result.stderr or result.stdout or "nvcc command failed").strip()
            return None, None, None, message

        match = re.search(r"release\s+(\d+)\.(\d+)", output)
        if not match:
            match = re.search(r"V(\d+)\.(\d+)", output)
        if not match:
            return None, None, None, "Unable to parse CUDA version from nvcc output"

        major = int(match.group(1))
        minor = int(match.group(2))
        return f"{major}.{minor}", major, minor, None

    def _get_distribution_version(self, names: List[str]) -> Optional[str]:
        for name in names:
            try:
                return metadata.version(name)
            except metadata.PackageNotFoundError:
                continue
        return None

    def _is_module_available(self, module_name: str) -> bool:
        importlib.invalidate_caches()
        if site.ENABLE_USER_SITE:
            try:
                user_site = site.getusersitepackages()
            except Exception:
                user_site = None
            if user_site and user_site not in sys.path:
                sys.path.append(user_site)
        if util.find_spec(module_name) is None:
            return False
        try:
            __import__(module_name)
        except Exception:
            return False
        return True

    def _is_site_writable(self) -> bool:
        try:
            site_packages = site.getsitepackages()
        except Exception:
            site_packages = []

        for path in site_packages:
            if path and os.path.isdir(path) and os.access(path, os.W_OK):
                return True
        return False

    def _truncate_output(self, output: Optional[str], limit: int = 4000) -> str:
        if not output:
            return ""
        if len(output) <= limit:
            return output
        return output[-limit:]

    def _get_recommended_pairs(self, cuda_major: int) -> List[Tuple[str, str]]:
        return CUPY_CUTENSOR_RECOMMENDED_BY_CUDA.get(cuda_major, [])

    def _get_libxc_requirement(self, cuda_major: int) -> str:
        cuda_suffix = f"{cuda_major}x"
        version = LIBXC_VERSION_BY_CUDA.get(cuda_major)
        if version:
            return f"gpu4pyscf-libxc-cuda{cuda_suffix}=={version}"
        return f"gpu4pyscf-libxc-cuda{cuda_suffix}"

    def _run_pip_install(
        self,
        packages: List[str],
        use_user_site: bool,
        extra_args: Optional[List[str]] = None
    ) -> Tuple[bool, str, str]:
        pip_command = [sys.executable, "-m", "pip", "install", "--no-cache-dir", "--prefer-binary"]
        if extra_args:
            pip_command.extend(extra_args)
        if use_user_site:
            pip_command.append("--user")
        pip_command.extend(packages)

        try:
            result = subprocess.run(
                pip_command,
                capture_output=True,
                text=True,
                timeout=PIP_INSTALL_TIMEOUT_SECONDS
            )
            stdout_tail = self._truncate_output(result.stdout)
            stderr_tail = self._truncate_output(result.stderr)
            return result.returncode == 0, stdout_tail, stderr_tail
        except subprocess.TimeoutExpired as exc:
            stdout_tail = self._truncate_output(exc.stdout)
            stderr_tail = self._truncate_output(exc.stderr)
            timeout_message = (
                f"pip install timed out after {PIP_INSTALL_TIMEOUT_SECONDS} seconds"
            )
            if stderr_tail:
                stderr_tail = f"{timeout_message}\n{stderr_tail}"
            else:
                stderr_tail = timeout_message
            return False, stdout_tail, stderr_tail

    def _build_pip_args(self, force_reinstall: bool, no_deps: bool = False) -> List[str]:
        args = ["--upgrade"]
        if force_reinstall:
            args.append("--force-reinstall")
        if no_deps:
            args.append("--no-deps")
        return args

    def _build_dependency_candidates(
        self,
        cuda_major: int,
        include_cutensor: bool
    ) -> List[List[str]]:
        cuda_suffix = f"{cuda_major}x"
        libxc_requirement = self._get_libxc_requirement(cuda_major)
        candidates: List[List[str]] = []

        for cupy_version, cutensor_version in self._get_recommended_pairs(cuda_major):
            deps = [
                f"cupy-cuda{cuda_suffix}=={cupy_version}",
                libxc_requirement,
            ]
            if include_cutensor:
                deps.append(f"cutensor-cu{cuda_major}=={cutensor_version}")
            candidates.append(deps)

        deps_latest = [
            f"cupy-cuda{cuda_suffix}",
            libxc_requirement,
        ]
        if include_cutensor:
            deps_latest.append(f"cutensor-cu{cuda_major}")
        candidates.append(deps_latest)

        return candidates

    def _install_dependency_first(
        self,
        gpu4pyscf_package: str,
        dependency_candidates: List[List[str]],
        use_user_site: bool,
        force_reinstall: bool
    ) -> Tuple[bool, List[str], str, str]:
        combined_stdout: List[str] = []
        combined_stderr: List[str] = []

        for dependency_packages in dependency_candidates:
            logger.info(f"Installing GPU4PySCF dependencies: {dependency_packages}")
            for use_no_deps in (False, True):
                attempt_stdout: List[str] = []
                attempt_stderr: List[str] = []
                dep_args = self._build_pip_args(force_reinstall, no_deps=use_no_deps)
                ok, stdout_tail, stderr_tail = self._run_pip_install(
                    dependency_packages,
                    use_user_site,
                    extra_args=dep_args
                )
                if stdout_tail:
                    attempt_stdout.append(stdout_tail)
                if stderr_tail:
                    attempt_stderr.append(stderr_tail)
                if not ok:
                    combined_stdout.extend(attempt_stdout)
                    combined_stderr.extend(attempt_stderr)
                    continue

                gpu_args = self._build_pip_args(force_reinstall, no_deps=True)
                ok_gpu, gpu_stdout, gpu_stderr = self._run_pip_install(
                    [gpu4pyscf_package],
                    use_user_site,
                    extra_args=gpu_args
                )
                if gpu_stdout:
                    attempt_stdout.append(gpu_stdout)
                if gpu_stderr:
                    attempt_stderr.append(gpu_stderr)
                if ok_gpu:
                    return True, dependency_packages + [gpu4pyscf_package], "\n\n".join(attempt_stdout), "\n\n".join(attempt_stderr)

                combined_stdout.extend(attempt_stdout)
                combined_stderr.extend(attempt_stderr)

        return False, [], "\n\n".join(combined_stdout), "\n\n".join(combined_stderr)

    def get_gpu4pyscf_status(self) -> Dict[str, Any]:
        is_linux = sys.platform.startswith("linux")
        cuda_version = None
        cuda_major = None
        cuda_minor = None
        detection_message = None

        if is_linux:
            cuda_version, cuda_major, cuda_minor, detection_message = self._detect_cuda_version()
        else:
            detection_message = "GPU4PySCF is supported on Linux only"

        cuda_detected = cuda_version is not None
        cuda_supported = cuda_major in CUDA_PACKAGE_MAP if cuda_major is not None else False

        recommended_gpu4pyscf = None
        recommended_cutensor = None
        if cuda_supported:
            recommended_gpu4pyscf, recommended_cutensor = CUDA_PACKAGE_MAP[cuda_major]
        elif cuda_detected:
            detection_message = detection_message or (
                f"Unsupported CUDA version {cuda_version}. Supported versions: 11.x, 12.x, 13.x"
            )

        gpu4pyscf_installed = self._is_module_available("gpu4pyscf")
        cutensor_installed = self._is_module_available("cutensor")

        gpu4pyscf_version = (
            self._get_distribution_version(
                ["gpu4pyscf-cuda13x", "gpu4pyscf-cuda12x", "gpu4pyscf-cuda11x", "gpu4pyscf"]
            )
            if gpu4pyscf_installed
            else None
        )
        cutensor_version = (
            self._get_distribution_version(
                ["cutensor-cu13", "cutensor-cu12", "cutensor-cu11", "cutensor"]
            )
            if cutensor_installed
            else None
        )

        return {
            "is_linux": is_linux,
            "cuda_detected": cuda_detected,
            "cuda_version": cuda_version,
            "cuda_major": cuda_major,
            "cuda_minor": cuda_minor,
            "cuda_supported": cuda_supported,
            "cuda_detection_message": detection_message,
            "recommended_gpu4pyscf_package": recommended_gpu4pyscf,
            "recommended_cutensor_package": recommended_cutensor,
            "gpu4pyscf_installed": gpu4pyscf_installed,
            "gpu4pyscf_version": gpu4pyscf_version,
            "cutensor_installed": cutensor_installed,
            "cutensor_version": cutensor_version,
        }

    def install_gpu4pyscf(
        self,
        include_cutensor: bool = True,
        force_reinstall: bool = False
    ) -> Dict[str, Any]:
        if not sys.platform.startswith("linux"):
            raise ValidationError("GPU4PySCF installation is supported on Linux only.")

        cuda_version, cuda_major, cuda_minor, detection_message = self._detect_cuda_version()
        if cuda_version is None or cuda_major is None:
            raise ValidationError(f"CUDA toolkit not detected: {detection_message}")

        if cuda_major not in CUDA_PACKAGE_MAP:
            raise ValidationError(
                f"Unsupported CUDA version {cuda_version}. Supported versions: 11.x, 12.x, 13.x"
            )

        gpu4pyscf_package, cutensor_package = CUDA_PACKAGE_MAP[cuda_major]
        packages = [gpu4pyscf_package]
        if include_cutensor:
            packages.append(cutensor_package)

        logger.info(f"Installing GPU4PySCF packages: {packages} (CUDA {cuda_version})")

        use_user_site = not self._is_site_writable()
        if use_user_site and not site.ENABLE_USER_SITE:
            raise ValidationError(
                "User site-packages is disabled; cannot install GPU4PySCF without a writable "
                "site-packages directory. Enable user site-packages or install into a writable "
                "environment."
            )
        dependency_candidates = self._build_dependency_candidates(
            cuda_major,
            include_cutensor
        )
        ok, installed_packages, stdout_tail, stderr_tail = self._install_dependency_first(
            gpu4pyscf_package,
            dependency_candidates,
            use_user_site,
            force_reinstall
        )

        if not ok:
            error_message = stderr_tail or stdout_tail or "pip install failed"
            if stdout_tail and stderr_tail:
                error_message = f"{error_message}\n\nFallback details:\n{stdout_tail}\n\n{stderr_tail}"
            raise ServiceError(f"pip install failed: {error_message}")

        status = self.get_gpu4pyscf_status()

        return {
            "status": status,
            "packages": installed_packages,
            "used_user_site": use_user_site,
            "pip_stdout": stdout_tail,
            "pip_stderr": stderr_tail,
        }
    
    def get_resource_status(self) -> Dict[str, Any]:
        """
        Get current system resource status including constraints and allocation.
        
        Returns:
            Dict containing system resource information
            
        Raises:
            ServiceError: If retrieval fails
        """
        try:
            logger.info("Getting system resource status")
            
            resource_manager = get_resource_manager()
            resource_summary = resource_manager.get_resource_summary()
            
            logger.info(f"Successfully retrieved system resource status")
            return resource_summary
            
        except Exception as e:
            logger.error(f"Failed to retrieve system resource status: {e}", exc_info=True)
            raise ServiceError(f'Failed to retrieve system resource status: {str(e)}')
    
    def get_system_diagnostics(self) -> Dict[str, Any]:
        """
        Get comprehensive system diagnostics for troubleshooting.
        
        Returns:
            Dict containing comprehensive diagnostic information
            
        Raises:
            ServiceError: If retrieval fails
        """
        try:
            logger.info("Getting comprehensive system diagnostics")
            
            # Load server config
            try:
                config_path = os.path.join(
                    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))),
                    'config',
                    'server-config.json'
                )
                if os.path.exists(config_path):
                    with open(config_path, 'r', encoding='utf-8') as f:
                        server_config = json.load(f)
                else:
                    server_config = {}
            except Exception:
                server_config = {}
            
            # Collect system information
            diagnostics = {
                'timestamp': datetime.now().isoformat(),
                'service_info': {
                    'service': 'pyscf-front-api',
                    'version': server_config.get('app_info', {}).get('version', 'unknown'),
                    'pid': os.getpid(),
                    'working_directory': os.getcwd()
                },
                'system_info': {
                    'cpu_count': multiprocessing.cpu_count(),
                    'platform': os.name,
                    'python_version': sys.version
                }
            }
            
            # Process manager diagnostics
            try:
                process_manager = get_process_manager()
                diagnostics['process_manager'] = {
                    'status': 'available',
                    'max_workers': process_manager.max_workers,
                    'max_parallel_instances': process_manager.max_parallel_instances,
                    'active_calculations': len(process_manager.active_futures),
                    'queued_calculations': len(process_manager.calculation_queue),
                    'is_shutdown': process_manager._shutdown,
                    'queue_status': process_manager.get_queue_status()
                }
            except Exception as pm_error:
                diagnostics['process_manager'] = {
                    'status': 'error',
                    'error': str(pm_error),
                    'error_type': type(pm_error).__name__
                }
            
            # Resource manager diagnostics
            try:
                resource_manager = get_resource_manager()
                diagnostics['resource_manager'] = resource_manager.get_diagnostics()
            except Exception as rm_error:
                diagnostics['resource_manager'] = {
                    'status': 'error',
                    'error': str(rm_error),
                    'error_type': type(rm_error).__name__
                }
            
            # File manager diagnostics
            try:
                # Load current settings to get calculations directory
                settings = get_current_settings()
                file_manager = CalculationRepository(base_dir=settings.calculations_directory)
                base_dir = file_manager.get_base_directory()
                
                diagnostics['file_manager'] = {
                    'status': 'available',
                    'base_directory': base_dir,
                    'base_directory_exists': os.path.exists(base_dir),
                    'base_directory_writable': os.access(base_dir, os.W_OK) if os.path.exists(base_dir) else False
                }
                
                # Count calculation directories
                try:
                    calculations = file_manager.list_calculations()
                    diagnostics['file_manager']['total_calculations'] = len(calculations)
                    diagnostics['file_manager']['calculation_statuses'] = {}
                    
                    # Count by status
                    status_counts = {}
                    for calc in calculations[:20]:  # Limit to first 20 for performance
                        try:
                            status = file_manager.read_calculation_status(os.path.join(base_dir, calc['id']))
                            status_counts[status] = status_counts.get(status, 0) + 1
                        except Exception:
                            status_counts['unknown'] = status_counts.get('unknown', 0) + 1
                    
                    diagnostics['file_manager']['calculation_statuses'] = status_counts
                except Exception as calc_error:
                    diagnostics['file_manager']['calculations_error'] = str(calc_error)
                    
            except Exception as fm_error:
                diagnostics['file_manager'] = {
                    'status': 'error',
                    'error': str(fm_error),
                    'error_type': type(fm_error).__name__
                }
            
            # Settings diagnostics
            try:
                settings = get_current_settings()
                diagnostics['settings'] = {
                    'status': 'available',
                    'settings': settings.model_dump()
                }
            except Exception as settings_error:
                diagnostics['settings'] = {
                    'status': 'error',
                    'error': str(settings_error),
                    'error_type': type(settings_error).__name__
                }
            
            logger.info("Successfully retrieved comprehensive system diagnostics")
            return diagnostics
            
        except Exception as e:
            logger.error(f"Failed to retrieve system diagnostics: {e}", exc_info=True)
            raise ServiceError(f'Failed to retrieve system diagnostics: {str(e)}')
    
    def get_process_manager_diagnostics(self) -> Dict[str, Any]:
        """
        Get detailed process manager diagnostics.
        
        Returns:
            Dict containing process manager diagnostic information
            
        Raises:
            ServiceError: If retrieval fails
        """
        try:
            logger.info("Getting process manager diagnostics")
            
            try:
                process_manager = get_process_manager()
                
                # Get detailed process manager state
                diagnostics = {
                    'timestamp': datetime.now().isoformat(),
                    'status': 'available',
                    'configuration': {
                        'max_workers': process_manager.max_workers,
                        'max_parallel_instances': process_manager.max_parallel_instances,
                        'is_shutdown': process_manager._shutdown
                    },
                    'current_state': {
                        'active_futures_count': len(process_manager.active_futures),
                        'active_calculation_ids': list(process_manager.active_futures.keys()),
                        'queued_calculations_count': len(process_manager.calculation_queue),
                        'completion_callbacks_count': len(process_manager.completion_callbacks)
                    },
                    'queue_details': [],
                    'resource_monitoring': {
                        'monitoring_active': process_manager._resource_monitor_thread is not None and process_manager._resource_monitor_thread.is_alive(),
                        'monitoring_interval': process_manager._resource_monitor_interval
                    }
                }
                
                # Get detailed queue information
                for i, queued_calc in enumerate(process_manager.calculation_queue[:10]):  # Limit to first 10
                    queue_item = {
                        'position': i + 1,
                        'calculation_id': queued_calc.calculation_id,
                        'created_at': queued_calc.created_at.isoformat(),
                        'waiting_reason': queued_calc.waiting_reason,
                        'calculation_method': queued_calc.parameters.get('calculation_method', 'unknown'),
                        'cpu_cores': queued_calc.parameters.get('cpu_cores', 'unknown'),
                        'memory_mb': queued_calc.parameters.get('memory_mb', 'unknown')
                    }
                    diagnostics['queue_details'].append(queue_item)
                
                # Executor status
                if process_manager.executor is not None:
                    diagnostics['executor'] = {
                        'available': True,
                        'type': type(process_manager.executor).__name__
                    }
                else:
                    diagnostics['executor'] = {
                        'available': False,
                        'error': 'ProcessPoolExecutor is None'
                    }
                    
            except Exception as pm_error:
                diagnostics = {
                    'timestamp': datetime.now().isoformat(),
                    'status': 'error',
                    'error': str(pm_error),
                    'error_type': type(pm_error).__name__
                }
            
            logger.info("Successfully retrieved process manager diagnostics")
            return diagnostics
            
        except Exception as e:
            logger.error(f"Failed to retrieve process manager diagnostics: {e}", exc_info=True)
            raise ServiceError(f'Failed to retrieve process manager diagnostics: {str(e)}')
    
    def get_resource_manager_diagnostics(self) -> Dict[str, Any]:
        """
        Get detailed resource manager diagnostics.
        
        Returns:
            Dict containing resource manager diagnostic information
            
        Raises:
            ServiceError: If retrieval fails
        """
        try:
            logger.info("Getting resource manager diagnostics")
            
            try:
                resource_manager = get_resource_manager()
                diagnostics = resource_manager.get_diagnostics()
                diagnostics['timestamp'] = datetime.now().isoformat()
            except Exception as rm_error:
                diagnostics = {
                    'timestamp': datetime.now().isoformat(),
                    'status': 'error',
                    'error': str(rm_error),
                    'error_type': type(rm_error).__name__
                }
            
            logger.info("Successfully retrieved resource manager diagnostics")
            return diagnostics
            
        except Exception as e:
            logger.error(f"Failed to retrieve resource manager diagnostics: {e}", exc_info=True)
            raise ServiceError(f'Failed to retrieve resource manager diagnostics: {str(e)}')
