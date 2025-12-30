"""
GPU management API endpoints.
Provides GPU status detection and gpu4pyscf installation helpers.
"""

import logging
from flask import Blueprint, jsonify, request
from flask_pydantic import validate

from services import get_gpu_service, ServiceError, ValidationError
from generated_models import GpuInstallRequest

logger = logging.getLogger(__name__)

# Create GPU blueprint
gpu_bp = Blueprint('gpu', __name__)


@gpu_bp.route('/api/gpu/status', methods=['GET'])
def get_gpu_status():
    """Get current GPU status and CUDA detection results."""
    try:
        gpu_service = get_gpu_service()
        status = gpu_service.get_status()
        return jsonify({'success': True, 'data': {'status': status}})
    except ServiceError as e:
        logger.error(f"Service error retrieving GPU status: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Failed to retrieve GPU status: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@gpu_bp.route('/api/gpu/install', methods=['POST'])
@validate()
def install_gpu(body: GpuInstallRequest):
    """Enqueue gpu4pyscf installation for detected CUDA toolkit."""
    try:
        gpu_service = get_gpu_service()
        job = gpu_service.enqueue_install_job(
            package=body.package, enable_gpu=body.enable_gpu
        )
        return jsonify({'success': True, 'data': {'job': job}}), 202
    except ValidationError as e:
        logger.warning(f"GPU install validation failed: {e.message}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except ServiceError as e:
        logger.error(f"Service error installing gpu4pyscf: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Failed to enqueue gpu4pyscf install: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@gpu_bp.route('/api/gpu/install/<job_id>', methods=['GET'])
def get_install_job(job_id: str):
    """Retrieve gpu4pyscf installation job status."""
    try:
        gpu_service = get_gpu_service()
        job = gpu_service.get_install_job(job_id)
        return jsonify({'success': True, 'data': {'job': job}})
    except ValidationError as e:
        logger.warning(f"GPU install job retrieval failed: {e.message}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except ServiceError as e:
        logger.error(f"Service error getting gpu install job: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Failed to fetch gpu install job: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@gpu_bp.route('/api/gpu/install/<job_id>/cancel', methods=['POST'])
def cancel_install_job(job_id: str):
    """Cancel a queued or running gpu4pyscf installation job."""
    try:
        gpu_service = get_gpu_service()
        job = gpu_service.cancel_install_job(job_id)
        return jsonify({'success': True, 'data': {'job': job}})
    except ValidationError as e:
        logger.warning(f"GPU install cancel validation failed: {e.message}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except ServiceError as e:
        logger.error(f"Service error cancelling gpu4pyscf job: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Failed to cancel gpu4pyscf job: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@gpu_bp.route('/api/gpu/install/jobs', methods=['GET'])
def list_install_jobs():
    """List recent gpu4pyscf installation jobs."""
    try:
        gpu_service = get_gpu_service()

        status_param = request.args.get('status')
        statuses = (
            [s.strip() for s in status_param.split(',') if s.strip()]
            if status_param
            else None
        )

        limit_param = request.args.get('limit')
        limit = 50
        if limit_param:
            try:
                limit = int(limit_param)
            except ValueError:
                return jsonify({'success': False, 'error': 'limit must be an integer'}), 400

        jobs = gpu_service.list_install_jobs(statuses=statuses, limit=limit)
        return jsonify({'success': True, 'data': {'jobs': jobs}})
    except ValidationError as e:
        logger.warning(f"GPU install job listing validation failed: {e.message}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except ServiceError as e:
        logger.error(f"Service error listing gpu4pyscf jobs: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Failed to list gpu4pyscf jobs: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500
