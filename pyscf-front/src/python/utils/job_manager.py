"""
Job Manager Utility
Manages calculation jobs and queue for the PySCF backend.
"""

import logging
from typing import Dict, Any, List, Optional
from enum import Enum
import threading
import queue
import time
from dataclasses import dataclass, field
from datetime import datetime

logger = logging.getLogger(__name__)


class JobStatus(Enum):
    """Enumeration of job statuses."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class Job:
    """Represents a calculation job."""
    id: str
    molecule_data: Dict[str, Any]
    calculation_params: Dict[str, Any]
    status: JobStatus = JobStatus.PENDING
    created_at: datetime = field(default_factory=datetime.now)
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    progress: float = 0.0


class JobManager:
    """Manages calculation jobs and execution queue."""
    
    def __init__(self, max_concurrent_jobs: int = 1):
        self.max_concurrent_jobs = max_concurrent_jobs
        self.jobs: Dict[str, Job] = {}
        self.job_queue = queue.Queue()
        self.running_jobs: Dict[str, Job] = {}
        self.worker_threads: List[threading.Thread] = []
        self.shutdown_event = threading.Event()
        
        # Start worker threads
        for i in range(max_concurrent_jobs):
            worker = threading.Thread(
                target=self._worker_thread,
                name=f"JobWorker-{i}",
                daemon=True
            )
            worker.start()
            self.worker_threads.append(worker)
        
        logger.info(f"JobManager initialized with {max_concurrent_jobs} worker threads")
    
    def submit_job(
        self,
        job_id: str,
        molecule_data: Dict[str, Any],
        calculation_params: Dict[str, Any]
    ) -> str:
        """
        Submit a new calculation job.
        
        Args:
            job_id: Unique identifier for the job
            molecule_data: Molecule structure data
            calculation_params: Calculation parameters
            
        Returns:
            Job ID
        """
        job = Job(
            id=job_id,
            molecule_data=molecule_data,
            calculation_params=calculation_params
        )
        
        self.jobs[job_id] = job
        self.job_queue.put(job)
        
        logger.info(f"Job {job_id} submitted to queue")
        return job_id
    
    def get_job_status(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Get status of a specific job."""
        if job_id not in self.jobs:
            return None
        
        job = self.jobs[job_id]
        return {
            'id': job.id,
            'status': job.status.value,
            'progress': job.progress,
            'created_at': job.created_at.isoformat(),
            'started_at': job.started_at.isoformat() if job.started_at else None,
            'completed_at': job.completed_at.isoformat() if job.completed_at else None,
            'result': job.result,
            'error': job.error
        }
    
    def get_status(self) -> Dict[str, Any]:
        """Get overall status of the job manager."""
        pending_count = sum(1 for job in self.jobs.values() if job.status == JobStatus.PENDING)
        running_count = len(self.running_jobs)
        completed_count = sum(1 for job in self.jobs.values() if job.status == JobStatus.COMPLETED)
        failed_count = sum(1 for job in self.jobs.values() if job.status == JobStatus.FAILED)
        
        return {
            'total_jobs': len(self.jobs),
            'pending': pending_count,
            'running': running_count,
            'completed': completed_count,
            'failed': failed_count,
            'queue_size': self.job_queue.qsize(),
            'max_concurrent_jobs': self.max_concurrent_jobs
        }
    
    def cancel_job(self, job_id: str) -> bool:
        """Cancel a job."""
        if job_id not in self.jobs:
            return False
        
        job = self.jobs[job_id]
        
        if job.status in [JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED]:
            return False
        
        job.status = JobStatus.CANCELLED
        job.completed_at = datetime.now()
        
        # Remove from running jobs if it's currently running
        if job_id in self.running_jobs:
            del self.running_jobs[job_id]
        
        logger.info(f"Job {job_id} cancelled")
        return True
    
    def _worker_thread(self):
        """Worker thread that processes jobs from the queue."""
        thread_name = threading.current_thread().name
        logger.info(f"Worker thread {thread_name} started")
        
        while not self.shutdown_event.is_set():
            try:
                # Get job from queue with timeout
                job = self.job_queue.get(timeout=1.0)
                
                if job.status == JobStatus.CANCELLED:
                    continue
                
                logger.info(f"Worker {thread_name} starting job {job.id}")
                self._execute_job(job)
                
            except queue.Empty:
                continue
            except Exception as e:
                logger.error(f"Worker {thread_name} encountered error: {e}")
    
    def _execute_job(self, job: Job):
        """Execute a single job."""
        try:
            # Mark job as running
            job.status = JobStatus.RUNNING
            job.started_at = datetime.now()
            self.running_jobs[job.id] = job
            
            # Import here to avoid circular imports
            from calculations.engine import CalculationEngine
            from utils.molecule_builder import MoleculeBuilder
            
            # Build molecule
            molecule_builder = MoleculeBuilder()
            mol = molecule_builder.build_from_data(job.molecule_data)
            
            # Perform calculation
            calculation_engine = CalculationEngine()
            result = calculation_engine.calculate(
                mol=mol,
                **job.calculation_params
            )
            
            # Mark job as completed
            job.status = JobStatus.COMPLETED
            job.completed_at = datetime.now()
            job.result = result
            job.progress = 1.0
            
            logger.info(f"Job {job.id} completed successfully")
            
        except Exception as e:
            # Mark job as failed
            job.status = JobStatus.FAILED
            job.completed_at = datetime.now()
            job.error = str(e)
            job.progress = 0.0
            
            logger.error(f"Job {job.id} failed: {e}")
            
        finally:
            # Remove from running jobs
            if job.id in self.running_jobs:
                del self.running_jobs[job.id]
    
    def shutdown(self):
        """Shutdown the job manager and all worker threads."""
        logger.info("Shutting down JobManager...")
        self.shutdown_event.set()
        
        # Wait for all worker threads to complete
        for worker in self.worker_threads:
            worker.join(timeout=5.0)
        
        logger.info("JobManager shutdown complete")