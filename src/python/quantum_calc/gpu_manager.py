"""GPU management utilities for gpu4pyscf integration on Linux."""

import logging
import platform
import queue
import re
import subprocess
import sys
import threading
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone, timedelta
from typing import Any, Dict, List, Optional, Tuple
from enum import Enum

try:
    from importlib import metadata
except ImportError:  # pragma: no cover - fallback for older Python
    import importlib_metadata as metadata  # type: ignore

from generated_models import GpuDevice, GpuStatus, Status as GpuReadyStatus

logger = logging.getLogger(__name__)


class GPUManagerError(Exception):
    """GPU manager level exception."""


class GpuInstallJobStatus(str, Enum):
    """Job status for gpu4pyscf installation queue."""

    queued = "queued"
    running = "running"
    succeeded = "succeeded"
    failed = "failed"
    canceled = "canceled"


@dataclass
class GPUInstallJobState:
    """State container for gpu4pyscf installation jobs."""

    job_id: str
    package: str
    enable_gpu: bool
    status: GpuInstallJobStatus = GpuInstallJobStatus.queued
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    logs: List[str] = field(default_factory=list)
    error: Optional[str] = None
    result_status: Optional[Dict[str, Any]] = None
    cancel_event: threading.Event = field(default_factory=threading.Event, repr=False)


class GPUManager:
    """Detects CUDA/GPU availability and manages gpu4pyscf installation."""

    SUPPORTED_METHODS = {"HF", "DFT", "TDDFT"}
    ALLOWED_PACKAGES = {
        "gpu4pyscf",
        "gpu4pyscf-cu11",
        "gpu4pyscf-cu12",
        "gpu4pyscf-cu13",
    }
    COMMAND_TIMEOUT_SEC = 5
    GPU_INSTALL_TIMEOUT_SEC = 900
    PROCESS_TERMINATE_TIMEOUT_SEC = 5
    JOB_RETENTION_SECONDS = 3600

    def __init__(self) -> None:
        self._last_status: Optional[GpuStatus] = None
        self._install_jobs: Dict[str, GPUInstallJobState] = {}
        self._job_queue: "queue.Queue[str]" = queue.Queue()
        self._worker_thread: Optional[threading.Thread] = None
        self._lock = threading.Lock()
        self._active_process: Optional[subprocess.Popen] = None
        self._active_job_id: Optional[str] = None

    def _start_worker(self):
        """Start background worker thread if not already running."""
        if self._worker_thread and self._worker_thread.is_alive():
            return

        self._worker_thread = threading.Thread(
            target=self._install_worker_loop, name="gpu-install-worker", daemon=True
        )
        self._worker_thread.start()

    def _parse_nvidia_smi_header(self, logs: List[str]) -> Tuple[Optional[str], Optional[str]]:
        """Parse CUDA/driver versions from `nvidia-smi` header."""
        try:
            output = subprocess.check_output(
                ["nvidia-smi"],
                text=True,
                timeout=self.COMMAND_TIMEOUT_SEC,
            )
            header = output.splitlines()[0] if output else ""
            cuda_match = re.search(r"CUDA Version:\s*([0-9.]+)", header)
            driver_match = re.search(r"Driver Version:\s*([0-9.]+)", header)
            return (
                cuda_match.group(1) if cuda_match else None,
                driver_match.group(1) if driver_match else None,
            )
        except subprocess.TimeoutExpired:
            logs.append(
                f"nvidia-smi timed out after {self.COMMAND_TIMEOUT_SEC}s"
            )
            return None, None
        except FileNotFoundError:
            logs.append("nvidia-smi not found on PATH")
            return None, None
        except subprocess.SubprocessError as e:
            logs.append(f"Failed to run nvidia-smi: {e}")
            return None, None

    def _query_gpu_list(self, logs: List[str]) -> Tuple[List[GpuDevice], bool]:
        """Return GPU list via nvidia-smi query."""
        gpus: List[GpuDevice] = []
        try:
            output = subprocess.check_output(
                [
                    "nvidia-smi",
                    "--query-gpu=name,memory.total",
                    "--format=csv,noheader,nounits",
                ],
                text=True,
                timeout=self.COMMAND_TIMEOUT_SEC,
            )
            for line in output.splitlines():
                parts = [p.strip() for p in line.split(",")]
                if not parts or not parts[0]:
                    continue
                memory_mb = None
                if len(parts) > 1:
                    try:
                        memory_mb = int(float(parts[1]))
                    except ValueError:
                        memory_mb = None
                gpus.append(GpuDevice(name=parts[0], memory_mb=memory_mb))
            return gpus, bool(gpus)
        except subprocess.TimeoutExpired:
            logs.append(
                f"nvidia-smi query timed out after {self.COMMAND_TIMEOUT_SEC}s"
            )
            return gpus, False
        except FileNotFoundError:
            logs.append("nvidia-smi not found on PATH")
            return gpus, False
        except subprocess.SubprocessError as e:
            logs.append(f"Failed to query GPU list via nvidia-smi: {e}")
            return gpus, False

    def _detect_cuda_version_with_nvcc(self, logs: List[str]) -> Optional[str]:
        """Fallback CUDA detection via nvcc."""
        try:
            output = subprocess.check_output(
                ["nvcc", "--version"],
                text=True,
                timeout=self.COMMAND_TIMEOUT_SEC,
            )
            match = re.search(r"release\s+([0-9.]+)", output)
            if match:
                return match.group(1)
        except subprocess.TimeoutExpired:
            logs.append(
                f"nvcc --version timed out after {self.COMMAND_TIMEOUT_SEC}s"
            )
        except FileNotFoundError:
            logs.append("nvcc not found on PATH")
        except subprocess.SubprocessError as e:
            logs.append(f"Failed to run nvcc --version: {e}")
        return None

    def _detect_gpu4pyscf(
        self, logs: List[str]
    ) -> Tuple[bool, Optional[str], Optional[str]]:
        """Check gpu4pyscf availability and version."""
        candidates = ["gpu4pyscf", "gpu4pyscf-cu13", "gpu4pyscf-cu12", "gpu4pyscf-cu11"]
        installed_package: Optional[str] = None
        installed_version: Optional[str] = None

        for name in candidates:
            try:
                installed_version = metadata.version(name)
                installed_package = name
                break
            except metadata.PackageNotFoundError:
                continue
            except Exception as e:  # pragma: no cover - unexpected metadata errors
                logs.append(f"Failed to inspect package metadata for {name}: {e}")

        if installed_package:
            try:
                import gpu4pyscf  # type: ignore

                installed_version = installed_version or getattr(
                    gpu4pyscf, "__version__", None
                )
                return True, installed_package, installed_version
            except Exception as e:
                logs.append(f"gpu4pyscf import failed: {e}")
                return False, installed_package, installed_version

        return False, None, None

    def _recommend_package(self, cuda_version: Optional[str]) -> Optional[str]:
        """Recommend gpu4pyscf package name based on CUDA version."""
        if not cuda_version:
            return None

        try:
            major = int(cuda_version.split(".")[0])
        except (ValueError, IndexError):
            return None

        if major >= 13:
            return "gpu4pyscf-cu13"
        if major >= 12:
            return "gpu4pyscf-cu12"
        if major == 11:
            return "gpu4pyscf-cu11"
        return None

    def get_status(self, force_refresh: bool = False) -> GpuStatus:
        """Return current GPU status (cached unless force_refresh)."""
        if self._last_status is not None and not force_refresh:
            return self._last_status

        logs: List[str] = []
        platform_name = platform.system()
        supported_platform = platform_name.lower() == "linux"

        cuda_version = None
        driver_version = None
        detected_gpus: List[GpuDevice] = []
        has_nvidia_gpu = False

        if supported_platform:
            cuda_version, driver_version = self._parse_nvidia_smi_header(logs)
            detected_gpus, has_nvidia_gpu = self._query_gpu_list(logs)
            if not cuda_version:
                cuda_version = self._detect_cuda_version_with_nvcc(logs)

        gpu4pyscf_installed, installed_package, gpu4pyscf_version = self._detect_gpu4pyscf(
            logs
        )
        recommended_package = self._recommend_package(cuda_version)

        # Determine readiness
        status_enum = GpuReadyStatus.error
        message = "GPU readiness unknown"
        if not supported_platform:
            status_enum = GpuReadyStatus.unsupported_platform
            message = "GPU acceleration is supported only on Linux."
        elif not has_nvidia_gpu:
            status_enum = GpuReadyStatus.missing_gpu
            message = "NVIDIA GPU was not detected. Install drivers and ensure nvidia-smi works."
        elif not cuda_version:
            status_enum = GpuReadyStatus.missing_cuda
            message = "CUDA version could not be detected. Verify NVIDIA driver and CUDA installation."
        elif not gpu4pyscf_installed:
            status_enum = GpuReadyStatus.not_installed
            message = "gpu4pyscf is not installed."
        else:
            status_enum = GpuReadyStatus.ready
            message = "GPU acceleration is ready."

        status = GpuStatus(
            platform=platform_name,
            gpu_supported_platform=supported_platform,
            has_nvidia_gpu=has_nvidia_gpu,
            cuda_version=cuda_version,
            driver_version=driver_version,
            detected_gpus=detected_gpus or None,
            gpu4pyscf_installed=gpu4pyscf_installed,
            gpu4pyscf_version=gpu4pyscf_version,
            installed_package=installed_package,
            recommended_package=recommended_package,
            status=status_enum,
            message=message,
            last_checked=datetime.now(timezone.utc),
            logs=logs or None,
        )

        self._last_status = status
        return status

    def _validate_package(self, package: str) -> str:
        """Ensure requested package is in the allowlist."""
        normalized = package.strip()
        if normalized not in self.ALLOWED_PACKAGES:
            allowed = ", ".join(sorted(self.ALLOWED_PACKAGES))
            raise GPUManagerError(f"Package not allowed. Allowed values: {allowed}")
        return normalized

    def enqueue_install_job(
        self, package: Optional[str] = None, enable_gpu: bool = True
    ) -> GPUInstallJobState:
        """
        Enqueue a gpu4pyscf installation job and return job metadata.

        This validates the environment up-front and then hands the work to a
        background worker thread so the API handler is not blocked.
        """
        status_snapshot = self.get_status(force_refresh=True)

        if not status_snapshot.gpu_supported_platform:
            raise GPUManagerError("GPU acceleration is only supported on Linux.")
        if not status_snapshot.has_nvidia_gpu:
            raise GPUManagerError(
                "NVIDIA GPU was not detected. Check your drivers and nvidia-smi."
            )
        if not status_snapshot.cuda_version:
            raise GPUManagerError(
                "CUDA was not detected. Verify your CUDA installation."
            )

        package_to_install = (
            package
            or status_snapshot.recommended_package
            or "gpu4pyscf"
        )
        package_to_install = self._validate_package(package_to_install)

        job_id = str(uuid.uuid4())
        job = GPUInstallJobState(
            job_id=job_id, package=package_to_install, enable_gpu=enable_gpu
        )
        job.logs.append(f"Queued gpu4pyscf install for package: {package_to_install}")

        with self._lock:
            self._install_jobs[job_id] = job
            self._job_queue.put(job_id)
            self._start_worker()

        return job

    def get_install_job(self, job_id: str) -> GPUInstallJobState:
        """Return install job state."""
        with self._lock:
            job = self._install_jobs.get(job_id)
        if job is None:
            raise GPUManagerError(f"Installation job not found: {job_id}")
        return job

    def _terminate_process(self, process: subprocess.Popen, job_logs: List[str], reason: str) -> None:
        """Terminate a subprocess safely with kill fallback."""
        job_logs.append(reason)
        try:
            process.terminate()
            process.wait(timeout=self.PROCESS_TERMINATE_TIMEOUT_SEC)
        except Exception as terminate_error:
            job_logs.append(f"Terminate failed, killing process: {terminate_error}")
            try:
                process.kill()
            except Exception as kill_error:
                job_logs.append(f"Kill failed: {kill_error}")

    def cancel_install_job(self, job_id: str) -> GPUInstallJobState:
        """Cancel a queued or running install job."""
        with self._lock:
            job = self._install_jobs.get(job_id)
            if job is None:
                raise GPUManagerError(f"Installation job not found: {job_id}")

            if job.status in {GpuInstallJobStatus.succeeded, GpuInstallJobStatus.failed, GpuInstallJobStatus.canceled}:
                return job

            job.cancel_event.set()
            if job.status == GpuInstallJobStatus.queued:
                job.status = GpuInstallJobStatus.canceled
                job.completed_at = datetime.now(timezone.utc)
                job.logs.append("Job canceled before start.")
                job.error = "Installation job was canceled before starting."
            elif (
                job.status == GpuInstallJobStatus.running
                and self._active_job_id == job_id
                and self._active_process is not None
            ):
                self._terminate_process(
                    self._active_process,
                    job.logs,
                    "Cancel requested. Terminating pip process...",
                )
        return job

    def list_install_jobs(
        self,
        statuses: Optional[List[GpuInstallJobStatus]] = None,
        limit: int = 50,
    ) -> List[GPUInstallJobState]:
        """
        List recent gpu4pyscf installation jobs (newest first).

        Args:
            statuses: Optional list of statuses to filter by.
            limit: Maximum number of jobs to return.
        """
        with self._lock:
            jobs = list(self._install_jobs.values())

        if statuses:
            status_set = set(statuses)
            jobs = [job for job in jobs if job.status in status_set]

        jobs.sort(key=lambda j: j.created_at, reverse=True)

        if limit > 0:
            jobs = jobs[:limit]

        return jobs

    def should_use_gpu_for_method(
        self, calculation_method: str
    ) -> Tuple[bool, str, GpuStatus]:
        """Return whether GPU should be used for a calculation method."""
        status = self.get_status(force_refresh=True)

        if not status.gpu_supported_platform:
            return False, "GPU acceleration is only supported on Linux.", status

        if status.status != GpuReadyStatus.ready:
            return (
                False,
                status.message or "GPU environment is not ready.",
                status,
            )

        if calculation_method not in self.SUPPORTED_METHODS:
            return (
                False,
                f"{calculation_method} is not GPU-enabled yet (supported: HF/DFT/TDDFT).",
                status,
            )

        return True, "GPU acceleration will be used.", status

    def try_create_gpu_hf(self, mol, spin: int):
        """
        Create gpu4pyscf HF object (RHF/UHF).

        Returns:
            Tuple of (hf_object or None, error_message or None)
        """
        try:
            from gpu4pyscf.scf import rhf, uhf  # type: ignore

            if spin == 0:
                return rhf.RHF(mol), None
            return uhf.UHF(mol), None
        except Exception as e:
            return None, str(e)

    def try_create_gpu_dft(self, mol, spin: int, xc: str):
        """
        Create gpu4pyscf DFT object (RKS/UKS).

        Returns:
            Tuple of (dft_object or None, error_message or None)
        """
        try:
            from gpu4pyscf.dft import rks, uks  # type: ignore

            if spin == 0:
                mf = rks.RKS(mol)
            else:
                mf = uks.UKS(mol)
            mf.xc = xc
            return mf, None
        except Exception as e:
            return None, str(e)

    def try_create_gpu_tddft(self, mf, tddft_method: str):
        """
        Create gpu4pyscf TDDFT/TDA object (RKS/UKS).

        Returns:
            Tuple of (tddft_object or None, error_message or None)
        """
        try:
            spin = getattr(mf.mol, "spin", 0)
            if tddft_method == "TDA":
                if spin == 0:
                    from gpu4pyscf.tdscf import rks as td_rks  # type: ignore
                    return td_rks.TDA(mf), None
                from gpu4pyscf.tdscf import uks as td_uks  # type: ignore
                return td_uks.TDA(mf), None

            if spin == 0:
                from gpu4pyscf.tdscf import rks as td_rks  # type: ignore
                return td_rks.TDDFT(mf), None
            from gpu4pyscf.tdscf import uks as td_uks  # type: ignore
            return td_uks.TDDFT(mf), None
        except Exception as e:
            return None, str(e)

    def _cleanup_completed_jobs(self) -> None:
        """Remove completed jobs older than retention window."""
        cutoff = datetime.now(timezone.utc) - timedelta(seconds=self.JOB_RETENTION_SECONDS)
        with self._lock:
            to_delete = [
                job_id
                for job_id, job in self._install_jobs.items()
                if job.completed_at and job.completed_at < cutoff
            ]
            for job_id in to_delete:
                self._install_jobs.pop(job_id, None)
        if to_delete:
            logger.info(
                f"Cleaned up {len(to_delete)} completed GPU install job(s) older than retention window"
            )

    def _install_worker_loop(self):
        """Process installation jobs from the queue."""
        while True:
            job_id = self._job_queue.get()
            with self._lock:
                job = self._install_jobs.get(job_id)
            if job is None:
                self._job_queue.task_done()
                continue

            if job.status == GpuInstallJobStatus.canceled:
                self._job_queue.task_done()
                continue

            with self._lock:
                job.status = GpuInstallJobStatus.running
                job.started_at = datetime.now(timezone.utc)
                self._active_job_id = job_id

            try:
                refreshed_status = self._perform_install_job(job)
                job.result_status = (
                    refreshed_status.model_dump(mode="json")
                    if hasattr(refreshed_status, "model_dump")
                    else None
                )
                if job.cancel_event.is_set():
                    job.error = job.error or "Installation was canceled."
                    job.logs.append("Job marked as canceled during install.")
                    job.status = GpuInstallJobStatus.canceled
                elif refreshed_status.status == GpuReadyStatus.ready:
                    job.logs.append("gpu4pyscf installation finished.")
                    job.status = GpuInstallJobStatus.succeeded
                else:
                    job.error = (
                        refreshed_status.message
                        or "GPU environment did not become ready."
                    )
                    job.logs.append(f"GPU not ready after install: {job.error}")
                    job.status = (
                        GpuInstallJobStatus.canceled
                        if job.cancel_event.is_set()
                        else GpuInstallJobStatus.failed
                    )
            except GPUManagerError as e:
                job.error = str(e)
                job.logs.append(str(e))
                job.status = (
                    GpuInstallJobStatus.canceled
                    if job.cancel_event.is_set()
                    else GpuInstallJobStatus.failed
                )
            except Exception as e:  # pragma: no cover - unexpected worker errors
                job.error = f"Failed to install gpu4pyscf: {e}"
                job.logs.append(job.error)
                job.status = GpuInstallJobStatus.failed
            finally:
                job.completed_at = datetime.now(timezone.utc)
                with self._lock:
                    self._active_process = None
                    self._active_job_id = None
                self._job_queue.task_done()
                self._cleanup_completed_jobs()

    def _perform_install_job(self, job: GPUInstallJobState) -> GpuStatus:
        """Execute pip install for a job (runs inside worker thread)."""
        status = self.get_status(force_refresh=True)
        if not status.gpu_supported_platform:
            raise GPUManagerError("GPU acceleration is only supported on Linux.")
        if not status.has_nvidia_gpu:
            raise GPUManagerError(
                "NVIDIA GPU was not detected. Check your drivers and nvidia-smi."
            )
        if job.cancel_event.is_set():
            raise GPUManagerError("Installation was canceled.")

        job.logs.append(
            f"Detected CUDA {status.cuda_version or 'unknown'}, recommended package: {status.recommended_package or 'gpu4pyscf'}"
        )
        job.logs.append(f"Installing package: {job.package}")

        timeout_timer: Optional[threading.Timer] = None
        try:
            process = subprocess.Popen(
                [sys.executable, "-m", "pip", "install", job.package],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            with self._lock:
                self._active_process = process

            def _on_timeout():
                job.cancel_event.set()
                with self._lock:
                    active_proc = self._active_process
                if active_proc and active_proc.poll() is None:
                    self._terminate_process(
                        active_proc,
                        job.logs,
                        f"Installation timed out after {self.GPU_INSTALL_TIMEOUT_SEC}s",
                    )

            timeout_timer = threading.Timer(self.GPU_INSTALL_TIMEOUT_SEC, _on_timeout)
            timeout_timer.daemon = True
            timeout_timer.start()
        except Exception as e:
            raise GPUManagerError(f"Failed to start pip: {e}") from e

        try:
            assert process.stdout is not None
            for line in process.stdout:
                cleaned = line.strip()
                if cleaned:
                    job.logs.append(cleaned)
                if job.cancel_event.is_set():
                    self._terminate_process(
                        process,
                        job.logs,
                        "Installation cancel requested; terminating pip process.",
                    )
                    raise GPUManagerError("Installation was canceled.")
            return_code = process.wait()
        finally:
            if timeout_timer:
                timeout_timer.cancel()
            with self._lock:
                self._active_process = None

        if job.cancel_event.is_set():
            raise GPUManagerError("Installation was canceled.")

        if return_code != 0:
            raise GPUManagerError("pip installation failed. See logs for details.")

        refreshed_status = self.get_status(force_refresh=True)
        combined_logs = (refreshed_status.logs or []) + job.logs
        refreshed_status.logs = combined_logs
        self._last_status = refreshed_status

        if job.enable_gpu and refreshed_status.status == GpuReadyStatus.ready:
            try:
                from quantum_calc import update_app_settings

                update_app_settings(
                    {
                        "gpu_acceleration_enabled": True,
                        "gpu_preferred_package": job.package,
                    }
                )
                refreshed_status.logs.append(
                    "GPU acceleration enabled in application settings."
                )
            except Exception as e:
                refreshed_status.logs.append(
                    f"Installed gpu4pyscf but failed to persist settings: {e}"
                )

        return refreshed_status

    def install_gpu4pyscf(
        self, package: Optional[str] = None, enable_gpu: bool = True
    ) -> GPUInstallJobState:
        """
        Backward-compatible synchronous entry point.

        Instead of blocking, it now enqueues and returns the job immediately.
        """
        return self.enqueue_install_job(package=package, enable_gpu=enable_gpu)


_gpu_manager: Optional[GPUManager] = None


def get_gpu_manager() -> GPUManager:
    """Get singleton GPUManager instance."""
    global _gpu_manager
    if _gpu_manager is None:
        _gpu_manager = GPUManager()
    return _gpu_manager
