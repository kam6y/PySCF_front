"""
System and debug API endpoints.
Handles system resource monitoring and diagnostic information.
"""

import logging
from datetime import datetime
from flask import Blueprint, jsonify

from services import get_system_service, ServiceError
from generated_models import SystemResourceSummary, SystemResourceInfo, ResourceConstraints, AllocatedResources

# Set up logging
logger = logging.getLogger(__name__)

# Create system blueprint
system_bp = Blueprint('system', __name__)


@system_bp.route('/api/system/resource-status', methods=['GET'])
def get_system_resource_status():
    """Get current system resource status including constraints and allocation."""
    try:
        system_service = get_system_service()
        
        # Call service layer
        resource_summary = system_service.get_resource_status()
        
        # Create response using Pydantic models
        system_info = SystemResourceInfo(
            total_cpu_cores=resource_summary['system_info']['total_cpu_cores'],
            total_memory_mb=resource_summary['system_info']['total_memory_mb'],
            available_memory_mb=resource_summary['system_info']['available_memory_mb'],
            cpu_usage_percent=resource_summary['system_info']['cpu_usage_percent'],
            memory_usage_percent=resource_summary['system_info']['memory_usage_percent'],
            timestamp=datetime.fromisoformat(resource_summary['system_info']['timestamp'].replace('Z', '+00:00'))
        )
        
        resource_constraints = ResourceConstraints(
            max_cpu_utilization_percent=resource_summary['resource_constraints']['max_cpu_utilization_percent'],
            max_memory_utilization_percent=resource_summary['resource_constraints']['max_memory_utilization_percent'],
            max_allowed_cpu_cores=resource_summary['resource_constraints']['max_allowed_cpu_cores'],
            max_allowed_memory_mb=resource_summary['resource_constraints']['max_allowed_memory_mb']
        )
        
        allocated_resources = AllocatedResources(
            total_allocated_cpu_cores=resource_summary['allocated_resources']['total_allocated_cpu_cores'],
            total_allocated_memory_mb=resource_summary['allocated_resources']['total_allocated_memory_mb'],
            available_cpu_cores=resource_summary['allocated_resources']['available_cpu_cores'],
            available_memory_mb=resource_summary['allocated_resources']['available_memory_mb'],
            active_calculations_count=resource_summary['allocated_resources']['active_calculations_count']
        )
        
        summary = SystemResourceSummary(
            system_info=system_info,
            resource_constraints=resource_constraints,
            allocated_resources=allocated_resources
        )
        
        return jsonify({
            'success': True,
            'data': summary.model_dump()
        })
        
    except ServiceError as e:
        logger.error(f"Service error retrieving system resource status: {e}")
        return jsonify({
            'success': False,
            'error': e.message
        }), e.status_code
    except Exception as e:
        logger.error(f"Failed to retrieve system resource status: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred.'
        }), 500


@system_bp.route('/api/debug/system-diagnostics', methods=['GET'])
def get_system_diagnostics():
    """Get comprehensive system diagnostics for troubleshooting."""
    try:
        system_service = get_system_service()
        
        # Call service layer
        diagnostics = system_service.get_system_diagnostics()
        
        return jsonify({
            'success': True,
            'data': diagnostics
        })
        
    except ServiceError as e:
        logger.error(f"Service error retrieving system diagnostics: {e}")
        return jsonify({
            'success': False,
            'error': e.message
        }), e.status_code
    except Exception as e:
        logger.error(f"Failed to retrieve system diagnostics: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred.'
        }), 500


@system_bp.route('/api/debug/process-manager-diagnostics', methods=['GET'])
def get_process_manager_diagnostics():
    """Get detailed process manager diagnostics."""
    try:
        system_service = get_system_service()
        
        # Call service layer
        diagnostics = system_service.get_process_manager_diagnostics()
        
        return jsonify({
            'success': True,
            'data': diagnostics
        })
        
    except ServiceError as e:
        logger.error(f"Service error retrieving process manager diagnostics: {e}")
        return jsonify({
            'success': False,
            'error': e.message
        }), e.status_code
    except Exception as e:
        logger.error(f"Failed to retrieve process manager diagnostics: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred.'
        }), 500


@system_bp.route('/api/debug/resource-manager-diagnostics', methods=['GET'])
def get_resource_manager_diagnostics():
    """Get detailed resource manager diagnostics."""
    try:
        system_service = get_system_service()
        
        # Call service layer
        diagnostics = system_service.get_resource_manager_diagnostics()
        
        return jsonify({
            'success': True,
            'data': diagnostics
        })
        
    except ServiceError as e:
        logger.error(f"Service error retrieving resource manager diagnostics: {e}")
        return jsonify({
            'success': False,
            'error': e.message
        }), e.status_code
    except Exception as e:
        logger.error(f"Failed to retrieve resource manager diagnostics: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred.'
        }), 500