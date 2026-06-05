"""
System resource monitoring service.

This service encapsulates system resource monitoring and diagnostic logic,
providing a unified interface for both API endpoints and AI agent tools.
"""

import logging
import os
import sys
import multiprocessing
import re
import importlib
import shutil
import site
import subprocess
from datetime import datetime
from typing import Dict, Any, Optional, Tuple, List
from importlib import metadata, util

from quantum_calc import (
    get_process_manager,
    get_current_settings,
    CalculationRepository,
    mask_settings,
)
from quantum_calc.resource_manager import get_resource_manager
from .exceptions import ServiceError, ValidationError
from config import get_server_config

logger = logging.getLogger(__name__)

PIP_INSTALL_TIMEOUT_SECONDS = 30 * 60
CUDA_DETECTION_TIMEOUT_SECONDS = 5

CUDA_PACKAGE_MAP = {
    11: ("gpu4pyscf-cuda11x", "cutensor-cu11"),
    12: ("gpu4pyscf-cuda12x", "cutensor-cu12"),
    13: ("gpu4pyscf-cuda13x", "cutensor-cu13"),
}

# Pinned versions for the top-level gpu4pyscf package per CUDA generation.
# These must be updated when upgrading to a new gpu4pyscf release.  Unlike
# cupy/cutensor, there is no version-fallback ladder — a wrong pin causes a
# total install failure.  Verify availability on PyPI before changing.
# The version is kept separate from CUDA_PACKAGE_MAP so that the status
# endpoint can display the package name without the ``==`` suffix.
#
# Verified against PyPI as of 2026-06-02:
#   gpu4pyscf-cuda11x: 0.6.1 .. 1.7.1
#   gpu4pyscf-cuda12x: 0.6.1 .. 1.7.1
#   gpu4pyscf-cuda13x: 1.4.3 .. 1.7.1
GPU4PYSCF_VERSION_BY_CUDA: dict[int, str] = {
    11: "1.7.1",
    12: "1.7.1",
    13: "1.7.1",
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

# Pinned fallback versions for cupy/cutensor when the recommended pair fails.
# These are used as a last-resort fallback and should be updated periodically.
CUPY_FALLBACK_VERSION = "13.4.1"
CUTENSOR_FALLBACK_VERSION = "2.2.0"

# Regex pattern to match URLs (http/https/ftp) so they can be preserved
# during path sanitization.
_URL_PATTERN = re.compile(r"https?://[^\s\"']+|ftp://[^\s\"']+", re.IGNORECASE)
# Regex pattern to detect and redact filesystem paths in pip output.
# Catches any absolute POSIX path (starting with /) that contains at least one
# path-separator-delimited component, covering standard directories as well as
# container/cloud mounts (/mnt, /data, /workspace, …) and virtualenv paths.
_PATH_SANITIZE_PATTERN = re.compile(
    r"(/(?:[A-Za-z0-9_.~-]+/)[^\s:\"']*)", re.IGNORECASE
)
# Matches Python traceback "File "/some/path.py"" lines and redacts the path.
_TRACEBACK_FILE_PATTERN = re.compile(r'(File\s+")[^"]*(/[^"]+)(")', re.IGNORECASE)
# Regex pattern to detect and redact environment variable references.
_ENV_SANITIZE_PATTERN = re.compile(
    r"(\b\w+_(?:PATH|HOME|DIR|KEY|TOKEN|SECRET)\s*=\s*)\S+", re.IGNORECASE
)


class SystemService:
    """Service for system resource monitoring and diagnostics."""

    def _detect_cuda_version(
        self,
    ) -> Tuple[Optional[str], Optional[int], Optional[int], Optional[str]]:
        if shutil.which("nvcc") is None:
            return None, None, None, "nvcc command not found"

        try:
            result = subprocess.run(
                ["nvcc", "--version"],
                capture_output=True,
                text=True,
                timeout=CUDA_DETECTION_TIMEOUT_SECONDS,
            )
            output = (result.stdout or "") + "\n" + (result.stderr or "")
        except subprocess.TimeoutExpired:
            return (
                None,
                None,
                None,
                f"nvcc --version timed out after {CUDA_DETECTION_TIMEOUT_SECONDS} seconds",
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

    @staticmethod
    def _sanitize_pip_output(text: str) -> str:
        """Remove local filesystem paths and env details from pip output.

        Redacts:
        - Absolute POSIX paths (any ``/dir/…`` with at least one component)
        - Python traceback ``File "/path/to/module.py"`` references
        - Environment variable assignments containing sensitive names

        URLs (http/https/ftp) are preserved to keep download progress and
        source information readable.
        """
        # 1. Temporarily replace URLs with numbered placeholders so the path
        #    regex does not redact URL paths.
        urls: list[str] = []

        def _stash_url(match: re.Match[str]) -> str:
            urls.append(match.group(0))
            return f"__URL_PLACEHOLDER_{len(urls) - 1}__"

        sanitized = _URL_PATTERN.sub(_stash_url, text)

        # 2. Redact Python traceback file references (more specific pattern).
        sanitized = _TRACEBACK_FILE_PATTERN.sub(r"\1<redacted-path>\3", sanitized)
        # 3. Redact remaining absolute POSIX paths.
        sanitized = _PATH_SANITIZE_PATTERN.sub("<redacted-path>", sanitized)
        # 4. Redact env-var assignments with sensitive-looking names.
        sanitized = _ENV_SANITIZE_PATTERN.sub(r"\1<redacted>", sanitized)

        # 5. Restore URLs.
        for i, url in enumerate(urls):
            sanitized = sanitized.replace(f"__URL_PLACEHOLDER_{i}__", url)

        return sanitized

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
        extra_args: Optional[List[str]] = None,
    ) -> Tuple[bool, str, str]:
        # NOTE: Residual supply-chain risk — transitive dependencies are NOT
        # hash-pinned.  Full ``--require-hashes`` is infeasible here because
        # the set of transitive wheels varies across CUDA variants and
        # platform ABIs, making a static hash-lock file impractical for a
        # runtime install.  Mitigations in place:
        #   1. Top-level packages are version-pinned (==).
        #   2. confirm_install gate prevents accidental installs.
        #   3. Loopback-address check restricts callers to the local machine.
        #   4. --no-input prevents pip from prompting (hangs).
        #   5. --disable-pip-version-check avoids leaking version info.
        #   6. No shell=True — argument list used throughout.
        pip_command = [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-cache-dir",
            "--prefer-binary",
            "--no-input",
            "--disable-pip-version-check",
        ]
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
                timeout=PIP_INSTALL_TIMEOUT_SECONDS,
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

    def _build_pip_args(
        self, force_reinstall: bool, no_deps: bool = False
    ) -> List[str]:
        args = ["--upgrade"]
        if force_reinstall:
            args.append("--force-reinstall")
        if no_deps:
            args.append("--no-deps")
        return args

    def _build_dependency_candidates(
        self, cuda_major: int, include_cutensor: bool
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

        # Pinned fallback: use known-good versions instead of unpinned latest.
        deps_fallback = [
            f"cupy-cuda{cuda_suffix}=={CUPY_FALLBACK_VERSION}",
            libxc_requirement,
        ]
        if include_cutensor:
            deps_fallback.append(
                f"cutensor-cu{cuda_major}=={CUTENSOR_FALLBACK_VERSION}"
            )
        candidates.append(deps_fallback)

        return candidates

    def _install_dependency_first(
        self,
        gpu4pyscf_package: str,
        dependency_candidates: List[List[str]],
        use_user_site: bool,
        force_reinstall: bool,
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
                    dependency_packages, use_user_site, extra_args=dep_args
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
                    [gpu4pyscf_package], use_user_site, extra_args=gpu_args
                )
                if gpu_stdout:
                    attempt_stdout.append(gpu_stdout)
                if gpu_stderr:
                    attempt_stderr.append(gpu_stderr)
                if ok_gpu:
                    return (
                        True,
                        dependency_packages + [gpu4pyscf_package],
                        "\n\n".join(attempt_stdout),
                        "\n\n".join(attempt_stderr),
                    )

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
            (
                cuda_version,
                cuda_major,
                cuda_minor,
                detection_message,
            ) = self._detect_cuda_version()
        else:
            detection_message = "GPU4PySCF is supported on Linux only"

        cuda_detected = cuda_version is not None
        cuda_supported = (
            cuda_major in CUDA_PACKAGE_MAP if cuda_major is not None else False
        )

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
                [
                    "gpu4pyscf-cuda13x",
                    "gpu4pyscf-cuda12x",
                    "gpu4pyscf-cuda11x",
                    "gpu4pyscf",
                ]
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

    @staticmethod
    def _build_gpu4pyscf_install_spec(package_name: str, cuda_major: int) -> str:
        """Build a version-pinned install spec for the gpu4pyscf package.

        Returns ``package_name==version`` when a pinned version exists for the
        given CUDA generation, or bare ``package_name`` as an unpinned fallback.
        """
        version = GPU4PYSCF_VERSION_BY_CUDA.get(cuda_major)
        if version is not None:
            return f"{package_name}=={version}"
        return package_name

    @staticmethod
    def _runtime_install_allowed() -> bool:
        """Check whether runtime package installation is permitted.

        Only known permissive environments (``development``, ``test``) allow
        runtime installs unconditionally.  All other values — including unset,
        empty, typos, and ``production`` — require the explicit opt-in
        ``PYSCF_ALLOW_RUNTIME_INSTALL=1``.  This mirrors the fail-closed
        pattern used by ``verify_auth_token`` in ``app.py``.
        """
        env = os.getenv("PYSCF_ENV", "").lower()
        if env in {"development", "test"}:
            return True
        raw = os.getenv("PYSCF_ALLOW_RUNTIME_INSTALL", "0")
        if raw not in {"0", "1", ""}:
            logger.warning(
                "Unrecognized PYSCF_ALLOW_RUNTIME_INSTALL=%r; only '1' enables runtime install",
                raw,
            )
        return raw == "1"

    def install_gpu4pyscf(
        self,
        include_cutensor: bool = True,
        force_reinstall: bool = False,
        confirm_install: bool = False,
    ) -> Dict[str, Any]:
        # In packaged builds, deny runtime installs unless explicitly opted in.
        if not self._runtime_install_allowed():
            raise ValidationError(
                "Runtime package installation is disabled in packaged builds. "
                "Set the PYSCF_ALLOW_RUNTIME_INSTALL=1 environment variable to "
                "enable it."
            )

        # Require explicit confirmation before mutating the Python environment.
        if not confirm_install:
            raise ValidationError(
                "GPU4PySCF installation mutates the Python environment. "
                "Set confirm_install=true to proceed."
            )

        if not sys.platform.startswith("linux"):
            raise ValidationError("GPU4PySCF installation is supported on Linux only.")

        (
            cuda_version,
            cuda_major,
            cuda_minor,
            detection_message,
        ) = self._detect_cuda_version()
        if cuda_version is None or cuda_major is None:
            raise ValidationError(f"CUDA toolkit not detected: {detection_message}")

        if cuda_major not in CUDA_PACKAGE_MAP:
            raise ValidationError(
                f"Unsupported CUDA version {cuda_version}. Supported versions: 11.x, 12.x, 13.x"
            )

        gpu4pyscf_package = CUDA_PACKAGE_MAP[cuda_major][0]

        # Pin the top-level gpu4pyscf package to a known-good version so that
        # ``pip install`` never silently pulls an unvetted release.
        gpu4pyscf_install_spec = self._build_gpu4pyscf_install_spec(
            gpu4pyscf_package, cuda_major
        )
        if "==" not in gpu4pyscf_install_spec:
            # Defensive fallback — if the CUDA generation is not in the map,
            # install unpinned (better than refusing entirely).
            logger.warning(
                "No pinned gpu4pyscf version for CUDA %d; installing unpinned",
                cuda_major,
            )

        logger.info(
            "Installing GPU4PySCF package: %s (CUDA %s, include_cutensor=%s)",
            gpu4pyscf_install_spec,
            cuda_version,
            include_cutensor,
        )

        use_user_site = not self._is_site_writable()
        if use_user_site and not site.ENABLE_USER_SITE:
            raise ValidationError(
                "User site-packages is disabled; cannot install GPU4PySCF without a writable "
                "site-packages directory. Enable user site-packages or install into a writable "
                "environment."
            )
        dependency_candidates = self._build_dependency_candidates(
            cuda_major, include_cutensor
        )
        (
            ok,
            installed_packages,
            stdout_tail,
            stderr_tail,
        ) = self._install_dependency_first(
            gpu4pyscf_install_spec,
            dependency_candidates,
            use_user_site,
            force_reinstall,
        )

        if not ok:
            # Sanitize error details before surfacing them
            safe_stderr = self._sanitize_pip_output(stderr_tail) if stderr_tail else ""
            safe_stdout = self._sanitize_pip_output(stdout_tail) if stdout_tail else ""
            error_message = safe_stderr or safe_stdout or "pip install failed"
            if safe_stdout and safe_stderr:
                error_message = f"{error_message}\n\nFallback details:\n{safe_stdout}\n\n{safe_stderr}"
            raise ServiceError(f"pip install failed: {error_message}")

        status = self.get_gpu4pyscf_status()

        # Sanitize pip output before returning to the client
        return {
            "status": status,
            "packages": installed_packages,
            "used_user_site": use_user_site,
            "pip_stdout": self._sanitize_pip_output(stdout_tail),
            "pip_stderr": self._sanitize_pip_output(stderr_tail),
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

            logger.info("Successfully retrieved system resource status")
            return resource_summary

        except Exception as e:
            logger.error(
                f"Failed to retrieve system resource status: {e}", exc_info=True
            )
            raise ServiceError(f"Failed to retrieve system resource status: {str(e)}")

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

            try:
                current_settings = get_current_settings()
            except Exception:
                current_settings = None

            try:
                _srv_cfg = get_server_config()
                app_version = _srv_cfg.get("app_info.version", "unknown")
            except Exception:
                app_version = "unknown"

            # Collect system information
            diagnostics = {
                "timestamp": datetime.now().isoformat(),
                "service_info": {
                    "service": "pyscf-front-api",
                    "version": app_version,
                    "pid": os.getpid(),
                    # Only expose the directory basename to avoid leaking the
                    # full filesystem path to API callers.
                    "working_directory": os.path.basename(os.getcwd()),
                },
                "system_info": {
                    "cpu_count": multiprocessing.cpu_count(),
                    "platform": os.name,
                    "python_version": sys.version,
                },
            }

            # Process manager diagnostics
            try:
                process_manager = get_process_manager()
                pm_diag = process_manager.get_diagnostics()
                diagnostics["process_manager"] = {
                    "status": "available",
                    "max_workers": pm_diag["max_workers"],
                    "max_parallel_instances": pm_diag["max_parallel_instances"],
                    "active_calculations": pm_diag["active_futures_count"],
                    "queued_calculations": pm_diag["queued_calculations_count"],
                    "is_shutdown": pm_diag["is_shutdown"],
                    "queue_status": process_manager.get_queue_status(),
                }
            except Exception as pm_error:
                diagnostics["process_manager"] = {
                    "status": "error",
                    "error": str(pm_error),
                    "error_type": type(pm_error).__name__,
                }

            # Resource manager diagnostics
            try:
                resource_manager = get_resource_manager()
                diagnostics["resource_manager"] = resource_manager.get_diagnostics()
            except Exception as rm_error:
                diagnostics["resource_manager"] = {
                    "status": "error",
                    "error": str(rm_error),
                    "error_type": type(rm_error).__name__,
                }

            # File manager diagnostics
            try:
                # Load current settings to get calculations directory
                if current_settings is None:
                    raise RuntimeError("Current settings unavailable")

                file_manager = CalculationRepository(
                    base_dir=current_settings.calculations_directory
                )
                base_dir = file_manager.get_base_directory()

                diagnostics["file_manager"] = {
                    "status": "available",
                    # Only expose the directory basename to avoid leaking the
                    # full filesystem path to API callers.
                    "base_directory": os.path.basename(base_dir),
                    "base_directory_exists": os.path.exists(base_dir),
                    "base_directory_writable": os.access(base_dir, os.W_OK)
                    if os.path.exists(base_dir)
                    else False,
                }

                # Count calculation directories
                try:
                    calculations = file_manager.list_calculations()
                    diagnostics["file_manager"]["total_calculations"] = len(
                        calculations
                    )
                    diagnostics["file_manager"]["calculation_statuses"] = {}

                    # Count by status
                    status_counts = {}
                    for calc in calculations[:20]:  # Limit to first 20 for performance
                        try:
                            status = file_manager.read_calculation_status(
                                os.path.join(base_dir, calc["id"])
                            )
                            status_counts[status] = status_counts.get(status, 0) + 1
                        except Exception:
                            status_counts["unknown"] = (
                                status_counts.get("unknown", 0) + 1
                            )

                    diagnostics["file_manager"]["calculation_statuses"] = status_counts
                except Exception as calc_error:
                    diagnostics["file_manager"]["calculations_error"] = str(calc_error)

            except Exception as fm_error:
                diagnostics["file_manager"] = {
                    "status": "error",
                    "error": str(fm_error),
                    "error_type": type(fm_error).__name__,
                }

            # Settings diagnostics
            try:
                if current_settings is None:
                    raise RuntimeError("Current settings unavailable")

                diagnostics["settings"] = {
                    "status": "available",
                    "settings": mask_settings(current_settings),
                }
            except Exception as settings_error:
                diagnostics["settings"] = {
                    "status": "error",
                    "error": str(settings_error),
                    "error_type": type(settings_error).__name__,
                }

            logger.info("Successfully retrieved comprehensive system diagnostics")
            return diagnostics

        except Exception as e:
            logger.error(f"Failed to retrieve system diagnostics: {e}", exc_info=True)
            raise ServiceError(f"Failed to retrieve system diagnostics: {str(e)}")

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
                pm_diag = process_manager.get_diagnostics()

                # Get detailed process manager state
                diagnostics = {
                    "timestamp": datetime.now().isoformat(),
                    "status": "available",
                    "configuration": {
                        "max_workers": pm_diag["max_workers"],
                        "max_parallel_instances": pm_diag["max_parallel_instances"],
                        "is_shutdown": pm_diag["is_shutdown"],
                    },
                    "current_state": {
                        "active_futures_count": pm_diag["active_futures_count"],
                        "active_calculation_ids": pm_diag["active_calculation_ids"],
                        "queued_calculations_count": pm_diag[
                            "queued_calculations_count"
                        ],
                    },
                    "queue_details": [],
                    "resource_monitoring": pm_diag["resource_monitoring"],
                }

                # Get detailed queue information
                for i, queued_calc in enumerate(
                    process_manager.calculation_queue[:10]
                ):  # Limit to first 10
                    queue_item = {
                        "position": i + 1,
                        "calculation_id": queued_calc.calculation_id,
                        "created_at": queued_calc.created_at.isoformat(),
                        "waiting_reason": queued_calc.waiting_reason,
                        "calculation_method": queued_calc.parameters.get(
                            "calculation_method", "unknown"
                        ),
                        "cpu_cores": queued_calc.parameters.get("cpu_cores", "unknown"),
                        "memory_mb": queued_calc.parameters.get("memory_mb", "unknown"),
                    }
                    diagnostics["queue_details"].append(queue_item)

                # Executor status
                diagnostics["executor"] = (
                    {
                        "available": pm_diag["executor_available"],
                        "type": pm_diag["executor_type"],
                    }
                    if pm_diag["executor_available"]
                    else {"available": False, "error": "ProcessPoolExecutor is None"}
                )

            except Exception as pm_error:
                diagnostics = {
                    "timestamp": datetime.now().isoformat(),
                    "status": "error",
                    "error": str(pm_error),
                    "error_type": type(pm_error).__name__,
                }

            logger.info("Successfully retrieved process manager diagnostics")
            return diagnostics

        except Exception as e:
            logger.error(
                f"Failed to retrieve process manager diagnostics: {e}", exc_info=True
            )
            raise ServiceError(
                f"Failed to retrieve process manager diagnostics: {str(e)}"
            )

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
                diagnostics["timestamp"] = datetime.now().isoformat()
            except Exception as rm_error:
                diagnostics = {
                    "timestamp": datetime.now().isoformat(),
                    "status": "error",
                    "error": str(rm_error),
                    "error_type": type(rm_error).__name__,
                }

            logger.info("Successfully retrieved resource manager diagnostics")
            return diagnostics

        except Exception as e:
            logger.error(
                f"Failed to retrieve resource manager diagnostics: {e}", exc_info=True
            )
            raise ServiceError(
                f"Failed to retrieve resource manager diagnostics: {str(e)}"
            )
