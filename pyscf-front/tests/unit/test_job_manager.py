"""
Unit tests for JobManager class
Tests job queue management and execution
"""

import unittest
import sys
import os
import time
import threading
from typing import Dict, Any

# Add src/python to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

from utils.job_manager import JobManager, JobStatus, Job
from utils.molecule_builder import MoleculeBuilder


class TestJobManager(unittest.TestCase):
    """Test cases for JobManager class"""
    
    def setUp(self):
        """Set up test fixtures"""
        # Use single worker to make tests predictable
        self.job_manager = JobManager(max_concurrent_jobs=1)
        self.builder = MoleculeBuilder()
        
        # Test molecule data
        self.test_molecule = self.builder.get_test_molecules()['hydrogen']
        self.test_calc_params = {
            'method': 'HF',
            'basis': 'STO-3G'
        }
    
    def tearDown(self):
        """Clean up after tests"""
        self.job_manager.shutdown()
    
    def test_job_submission(self):
        """Test basic job submission"""
        job_id = self.job_manager.submit_job(
            job_id='test_job_1',
            molecule_data=self.test_molecule,
            calculation_params=self.test_calc_params
        )
        
        self.assertEqual(job_id, 'test_job_1')
        
        # Check job exists in manager
        job_status = self.job_manager.get_job_status('test_job_1')
        self.assertIsNotNone(job_status)
        self.assertEqual(job_status['id'], 'test_job_1')
        self.assertEqual(job_status['status'], JobStatus.PENDING.value)
    
    def test_job_execution(self):
        """Test job execution and completion"""
        job_id = self.job_manager.submit_job(
            job_id='test_job_2',
            molecule_data=self.test_molecule,
            calculation_params=self.test_calc_params
        )
        
        # Wait for job to complete (should be quick for H2/STO-3G)
        max_wait = 10  # seconds
        start_time = time.time()
        
        while time.time() - start_time < max_wait:
            status = self.job_manager.get_job_status(job_id)
            if status['status'] in [JobStatus.COMPLETED.value, JobStatus.FAILED.value]:
                break
            time.sleep(0.1)
        
        # Check final status
        final_status = self.job_manager.get_job_status(job_id)
        self.assertEqual(final_status['status'], JobStatus.COMPLETED.value)
        
        # Should have results
        self.assertIsNotNone(final_status['result'])
        self.assertIn('energy', final_status['result'])
        self.assertIn('method', final_status['result'])
    
    def test_multiple_jobs(self):
        """Test submitting multiple jobs"""
        job_ids = []
        
        for i in range(3):
            job_id = f'test_job_multi_{i}'
            self.job_manager.submit_job(
                job_id=job_id,
                molecule_data=self.test_molecule,
                calculation_params=self.test_calc_params
            )
            job_ids.append(job_id)
        
        # Wait for all jobs to complete
        max_wait = 30  # seconds
        start_time = time.time()
        
        while time.time() - start_time < max_wait:
            statuses = [self.job_manager.get_job_status(jid) for jid in job_ids]
            completed = [s['status'] == JobStatus.COMPLETED.value for s in statuses]
            
            if all(completed):
                break
            time.sleep(0.2)
        
        # Check all jobs completed
        for job_id in job_ids:
            status = self.job_manager.get_job_status(job_id)
            self.assertEqual(status['status'], JobStatus.COMPLETED.value)
            self.assertIsNotNone(status['result'])
    
    def test_job_cancellation(self):
        """Test job cancellation"""
        job_id = self.job_manager.submit_job(
            job_id='test_job_cancel',
            molecule_data=self.test_molecule,
            calculation_params=self.test_calc_params
        )
        
        # Cancel immediately
        success = self.job_manager.cancel_job(job_id)
        self.assertTrue(success)
        
        # Check status
        status = self.job_manager.get_job_status(job_id)
        self.assertEqual(status['status'], JobStatus.CANCELLED.value)
    
    def test_manager_status(self):
        """Test getting overall manager status"""
        # Submit some jobs
        for i in range(2):
            self.job_manager.submit_job(
                job_id=f'status_test_{i}',
                molecule_data=self.test_molecule,
                calculation_params=self.test_calc_params
            )
        
        status = self.job_manager.get_status()
        
        # Check status structure
        self.assertIn('total_jobs', status)
        self.assertIn('pending', status)
        self.assertIn('running', status)
        self.assertIn('completed', status)
        self.assertIn('failed', status)
        self.assertIn('queue_size', status)
        self.assertIn('max_concurrent_jobs', status)
        
        # Should have jobs
        self.assertGreaterEqual(status['total_jobs'], 2)
        self.assertEqual(status['max_concurrent_jobs'], 1)
    
    def test_nonexistent_job(self):
        """Test querying nonexistent job"""
        status = self.job_manager.get_job_status('nonexistent_job')
        self.assertIsNone(status)
    
    def test_cancel_nonexistent_job(self):
        """Test cancelling nonexistent job"""
        success = self.job_manager.cancel_job('nonexistent_job')
        self.assertFalse(success)
    
    def test_cancel_completed_job(self):
        """Test cancelling already completed job"""
        job_id = self.job_manager.submit_job(
            job_id='test_complete_cancel',
            molecule_data=self.test_molecule,
            calculation_params=self.test_calc_params
        )
        
        # Wait for completion
        max_wait = 10
        start_time = time.time()
        
        while time.time() - start_time < max_wait:
            status = self.job_manager.get_job_status(job_id)
            if status['status'] == JobStatus.COMPLETED.value:
                break
            time.sleep(0.1)
        
        # Try to cancel completed job
        success = self.job_manager.cancel_job(job_id)
        self.assertFalse(success)
    
    def test_job_timestamps(self):
        """Test job timestamp tracking"""
        job_id = self.job_manager.submit_job(
            job_id='test_timestamps',
            molecule_data=self.test_molecule,
            calculation_params=self.test_calc_params
        )
        
        # Check initial timestamp
        initial_status = self.job_manager.get_job_status(job_id)
        self.assertIsNotNone(initial_status['created_at'])
        self.assertIsNone(initial_status['started_at'])
        self.assertIsNone(initial_status['completed_at'])
        
        # Wait for completion
        max_wait = 10
        start_time = time.time()
        
        while time.time() - start_time < max_wait:
            status = self.job_manager.get_job_status(job_id)
            if status['status'] == JobStatus.COMPLETED.value:
                break
            time.sleep(0.1)
        
        # Check final timestamps
        final_status = self.job_manager.get_job_status(job_id)
        self.assertIsNotNone(final_status['started_at'])
        self.assertIsNotNone(final_status['completed_at'])
    
    def test_failed_job_handling(self):
        """Test handling of failed jobs"""
        # Submit job with invalid parameters to cause failure
        invalid_params = {
            'method': 'INVALID_METHOD',
            'basis': 'STO-3G'
        }
        
        job_id = self.job_manager.submit_job(
            job_id='test_failed_job',
            molecule_data=self.test_molecule,
            calculation_params=invalid_params
        )
        
        # Wait for job to fail
        max_wait = 10
        start_time = time.time()
        
        while time.time() - start_time < max_wait:
            status = self.job_manager.get_job_status(job_id)
            if status['status'] in [JobStatus.FAILED.value, JobStatus.COMPLETED.value]:
                break
            time.sleep(0.1)
        
        # Should have failed
        final_status = self.job_manager.get_job_status(job_id)
        self.assertEqual(final_status['status'], JobStatus.FAILED.value)
        self.assertIsNotNone(final_status['error'])
        self.assertIsNone(final_status['result'])


if __name__ == '__main__':
    unittest.main()