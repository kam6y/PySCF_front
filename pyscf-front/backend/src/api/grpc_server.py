"""gRPC Server Implementation for PySCF_Front"""

import logging
import asyncio
from typing import Optional, Dict, Any, Iterator
import time
import uuid

import grpc
from concurrent.futures import ThreadPoolExecutor

from core.pyscf_engine import engine, MoleculeData, AtomData
from core.database import DatabaseManager, Instance, Molecule, Calculation, Result
from generated import calculation_pb2_grpc, calculation_pb2

logger = logging.getLogger(__name__)


class CalculationServicer(calculation_pb2_grpc.CalculationServiceServicer):
    """gRPC service implementation for quantum chemistry calculations"""
    
    def __init__(self):
        self.active_calculations: Dict[str, Any] = {}
        self.executor = ThreadPoolExecutor(max_workers=4)
    
    def HealthCheck(self, request, context):
        """Health check endpoint"""
        try:
            pyscf_status = engine.health_check()
            
            response = calculation_pb2.HealthResponse()
            response.healthy = pyscf_status.get("available", False)
            response.status = "healthy" if response.healthy else "unhealthy"
            response.details = f"PySCF available: {pyscf_status.get('available', False)}"
            
            return response
            
        except Exception as e:
            logger.error(f"Health check failed: {e}")
            response = calculation_pb2.HealthResponse()
            response.healthy = False
            response.status = "error"
            response.details = str(e)
            return response
    
    def GetSystemInfo(self, request, context):
        """Get system information"""
        try:
            info = engine.get_system_info()
            
            response = calculation_pb2.SystemInfoResponse()
            response.success = True
            response.message = "System info retrieved successfully"
            
            # Fill system info
            system_info = response.system_info
            system_info.version = info["version"]
            system_info.python_version = info["python_version"]
            system_info.pyscf_version = info.get("pyscf_version", "Not available")
            system_info.gpu_available = info["gpu_available"]
            
            # Resources
            resources = system_info.resources
            resources.cpu_cores = info["resources"]["cpu_cores"]
            resources.total_memory_mb = info["resources"]["total_memory_mb"]
            
            # Available methods and basis sets
            system_info.available_methods.extend(info["available_methods"])
            system_info.available_basis_sets.extend(info["available_basis_sets"])
            
            return response
            
        except Exception as e:
            logger.error(f"System info failed: {e}")
            response = calculation_pb2.SystemInfoResponse()
            response.success = False
            response.message = str(e)
            return response
    
    def CreateMolecule(self, request, context):
        """Create a new molecule"""
        try:
            # Convert gRPC atoms to internal format
            atoms = []
            for grpc_atom in request.atoms:
                atoms.append(AtomData(
                    symbol=grpc_atom.symbol,
                    x=grpc_atom.x,
                    y=grpc_atom.y,
                    z=grpc_atom.z
                ))
            
            molecule_data = MoleculeData(
                name=request.name,
                formula=request.formula,
                atoms=atoms,
                charge=request.charge,
                multiplicity=request.multiplicity
            )
            
            # Create response
            response = calculation_pb2.MoleculeResponse()
            response.success = True
            response.message = "Molecule created successfully"
            
            # Fill molecule data
            molecule = response.molecule
            molecule.id = str(uuid.uuid4())
            molecule.name = molecule_data.name
            molecule.formula = molecule_data.formula
            molecule.molecular_weight = len(atoms) * 10.0  # Simple estimation
            molecule.charge = molecule_data.charge
            molecule.multiplicity = molecule_data.multiplicity
            
            # Add atoms to response
            for atom_data in atoms:
                atom = molecule.atoms.add()
                atom.symbol = atom_data.symbol
                atom.x = atom_data.x
                atom.y = atom_data.y
                atom.z = atom_data.z
            
            return response
            
        except Exception as e:
            logger.error(f"Create molecule failed: {e}")
            response = calculation_pb2.MoleculeResponse()
            response.success = False
            response.message = str(e)
            return response
    
    def CreateInstance(self, request, context):
        """Create a calculation instance"""
        try:
            # Create response
            response = calculation_pb2.InstanceResponse()
            response.success = True
            response.message = "Instance created successfully"
            
            # Fill instance data
            instance = response.instance
            instance.id = str(uuid.uuid4())
            instance.name = request.name
            instance.description = request.description
            
            return response
            
        except Exception as e:
            logger.error(f"Create instance failed: {e}")
            response = calculation_pb2.InstanceResponse()
            response.success = False
            response.message = str(e)
            return response
    
    def StartCalculation(self, request, context):
        """Start a calculation and stream progress updates"""
        calculation_id = str(uuid.uuid4())
        
        try:
            # Create a simple molecule for testing (in real implementation, get from database)
            molecule_data = MoleculeData(
                name="Test Molecule",
                formula="H2O",
                atoms=[
                    AtomData("O", 0.0, 0.0, 0.0),
                    AtomData("H", 0.757, 0.586, 0.0),
                    AtomData("H", -0.757, 0.586, 0.0)
                ]
            )
            
            # Start calculation in background
            calculation_generator = engine.run_calculation(
                molecule_data=molecule_data,
                method=request.method,
                basis_set=request.basis_set,
                max_iterations=100
            )
            
            # Stream progress updates
            for progress_data in calculation_generator:
                progress = calculation_pb2.CalculationProgress()
                progress.job_id = calculation_id
                progress.current_step = progress_data.get("step", "")
                progress.progress_percentage = float(progress_data.get("progress", 0))
                progress.message = progress_data.get("message", "")
                
                if progress_data.get("step") == "complete":
                    progress.status = calculation_pb2.CalculationStatus.CALCULATION_STATUS_COMPLETED
                else:
                    progress.status = calculation_pb2.CalculationStatus.CALCULATION_STATUS_RUNNING
                
                yield progress
            
            # Final completion message
            final_progress = calculation_pb2.CalculationProgress()
            final_progress.job_id = calculation_id
            final_progress.current_step = "completed"
            final_progress.progress_percentage = 100.0
            final_progress.message = "Calculation completed successfully"
            final_progress.status = calculation_pb2.CalculationStatus.CALCULATION_STATUS_COMPLETED
            
            yield final_progress
            
        except Exception as e:
            logger.error(f"Calculation failed: {e}")
            
            # Send error status
            error_progress = calculation_pb2.CalculationProgress()
            error_progress.job_id = calculation_id
            error_progress.current_step = "error"
            error_progress.progress_percentage = 0.0
            error_progress.message = f"Calculation failed: {str(e)}"
            error_progress.status = calculation_pb2.CalculationStatus.CALCULATION_STATUS_FAILED
            
            yield error_progress
    
    def GetResults(self, request, context):
        """Get calculation results"""
        try:
            # For now, return mock results
            response = calculation_pb2.ResultResponse()
            response.success = True
            response.message = "Results retrieved successfully"
            
            # Mock energy data
            energy_data = response.energy_data
            energy_data.total_energy = -76.026749  # Hartree
            energy_data.homo_energy = -0.491467
            energy_data.lumo_energy = 0.151735
            energy_data.homo_lumo_gap = 0.643202
            
            # Mock orbital data
            for i in range(5):  # 5 orbitals for water
                orbital = response.orbitals.add()
                orbital.index = i
                orbital.energy = -1.0 + i * 0.5
                orbital.occupancy = 2.0 if i < 2 else 0.0
                orbital.type = "homo" if i == 1 else ("lumo" if i == 2 else "occupied" if i < 2 else "virtual")
            
            return response
            
        except Exception as e:
            logger.error(f"Get results failed: {e}")
            response = calculation_pb2.ResultResponse()
            response.success = False
            response.message = str(e)
            return response
    
    def CancelCalculation(self, request, context):
        """Cancel a running calculation"""
        try:
            # Implementation for cancelling calculations
            response = calculation_pb2.CancelResponse()
            response.success = True
            response.message = "Calculation cancelled successfully"
            return response
            
        except Exception as e:
            logger.error(f"Cancel calculation failed: {e}")
            response = calculation_pb2.CancelResponse()
            response.success = False
            response.message = str(e)
            return response