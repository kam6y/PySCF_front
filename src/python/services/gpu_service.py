"""
GPU management service.

Provides GPU detection status and installation flow for gpu4pyscf.
"""

import logging
from typing import Any, Dict, List, Optional

from quantum_calc.gpu_manager import (
    GPUManagerError,
    GPUInstallJobState,
    GpuInstallJobStatus,
    get_gpu_manager,
)
from .exceptions import ServiceError, ValidationError

logger = logging.getLogger(__name__)


class GPUService:
    """Service for GPU detection and gpu4pyscf installation."""

    def __init__(self):
        self.gpu_manager = get_gpu_manager()

    def _serialize_job(self, job: GPUInstallJobState) -> Dict[str, Any]:
        """Convert internal job state to a JSON-safe dict."""
        return {
            "job_id": job.job_id,
            "package": job.package,
            "enable_gpu": job.enable_gpu,
            "status": job.status.value,
            "error": job.error,
            "logs": job.logs or None,
            "created_at": job.created_at.isoformat() if job.created_at else None,
            "started_at": job.started_at.isoformat() if job.started_at else None,
            "completed_at": job.completed_at.isoformat() if job.completed_at else None,
            "result_status": job.result_status,
        }

    def get_status(self) -> Dict[str, Any]:
        """Return current GPU status."""
        try:
            status = self.gpu_manager.get_status(force_refresh=True)
            return status.model_dump(mode="json")
        except Exception as e:
            logger.error(f"Failed to retrieve GPU status: {e}", exc_info=True)
            raise ServiceError(f"Failed to retrieve GPU status: {e}")

    def enqueue_install_job(
        self, package: Optional[str] = None, enable_gpu: bool = True
    ) -> Dict[str, Any]:
        """Enqueue gpu4pyscf installation and return job metadata."""
        try:
            job = self.gpu_manager.enqueue_install_job(
                package=package, enable_gpu=enable_gpu
            )
            return self._serialize_job(job)
        except GPUManagerError as e:
            logger.warning(f"GPU installation validation failed: {e}")
            raise ValidationError(str(e))
        except ValidationError:
            raise
        except Exception as e:
            logger.error(f"Failed to enqueue gpu4pyscf install: {e}", exc_info=True)
            raise ServiceError(f"Failed to enqueue gpu4pyscf install: {e}")

    def get_install_job(self, job_id: str) -> Dict[str, Any]:
        """Fetch install job state."""
        try:
            job = self.gpu_manager.get_install_job(job_id)
            return self._serialize_job(job)
        except GPUManagerError as e:
            raise ValidationError(str(e))
        except ValidationError:
            raise
        except Exception as e:
            logger.error(f"Failed to get gpu install job {job_id}: {e}", exc_info=True)
            raise ServiceError(f"Failed to get gpu install job: {e}")

    def cancel_install_job(self, job_id: str) -> Dict[str, Any]:
        """Cancel a queued/running install job."""
        try:
            job = self.gpu_manager.cancel_install_job(job_id)
            return self._serialize_job(job)
        except GPUManagerError as e:
            raise ValidationError(str(e))
        except ValidationError:
            raise
        except Exception as e:
            logger.error(f"Failed to cancel gpu install job {job_id}: {e}", exc_info=True)
            raise ServiceError(f"Failed to cancel gpu install job: {e}")

    def list_install_jobs(
        self,
        statuses: Optional[List[str]] = None,
        limit: int = 50,
    ) -> List[Dict[str, Any]]:
        """List recent install jobs, newest first."""
        try:
            parsed_statuses: Optional[List[GpuInstallJobStatus]] = None
            if statuses:
                parsed_statuses = []
                for status in statuses:
                    try:
                        parsed_statuses.append(GpuInstallJobStatus(status))
                    except ValueError:
                        raise ValidationError(f"Invalid job status: {status}")

            limit_clamped = max(1, min(limit, 200))
            jobs = self.gpu_manager.list_install_jobs(
                statuses=parsed_statuses, limit=limit_clamped
            )
            return [self._serialize_job(job) for job in jobs]
        except ValidationError:
            raise
        except Exception as e:
            logger.error(f"Failed to list gpu install jobs: {e}", exc_info=True)
            raise ServiceError(f"Failed to list gpu install jobs: {e}")
