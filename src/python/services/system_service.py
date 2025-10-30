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
from datetime import datetime
from typing import Dict, Any

from quantum_calc import get_process_manager, get_current_settings
from quantum_calc.resource_manager import get_resource_manager
from quantum_calc.file_manager import CalculationFileManager
from .exceptions import ServiceError

logger = logging.getLogger(__name__)


class SystemService:
    """Service for system resource monitoring and diagnostics."""
    
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
                file_manager = CalculationFileManager(base_dir=settings.calculations_directory)
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
