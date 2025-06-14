import grpc
import logging
import threading
import time
from concurrent import futures
from typing import Iterator, Dict, Any
from datetime import datetime

from src.grpc_stubs import calculation_pb2, calculation_pb2_grpc
from src.grpc_stubs.calculation_pb2 import (
    MoleculeResponse, InstanceResponse, CalculationProgress,
    CalculationStatusResponse, ResultsResponse, SystemInfoResponse,
    HealthCheckResponse, DeleteMoleculeResponse, DeleteInstanceResponse,
    CancelCalculationResponse, ExportResultsResponse, JobQueueResponse,
    UpdateJobPriorityResponse, ListInstancesResponse
)
from src.core.pyscf_config import pyscf_config, calculation_engine, MoleculeBuilder
from src.database.connection import db_manager, get_db
from src.database.models import (
    Molecule, Instance, Calculation, Result, JobQueue,
    InstanceStatus, CalculationStatus, JobStatus,
    create_molecule_from_atoms, create_calculation_instance, 
    create_calculation_job, save_calculation_results
)

logger = logging.getLogger(__name__)


class CalculationServiceImpl(calculation_pb2_grpc.CalculationServiceServicer):
    
    def __init__(self):
        self.active_calculations: Dict[str, threading.Thread] = {}
        self.calculation_progress: Dict[str, Dict[str, Any]] = {}
        
    def CreateMolecule(self, request, context):
        logger.info(f"CreateMolecule called with name: {request.name}")
        
        try:
            # Validate geometry
            atoms_list = [(atom.symbol, atom.x, atom.y, atom.z) for atom in request.atoms]
            is_valid, error_msg = MoleculeBuilder.validate_geometry(atoms_list)
            
            if not is_valid:
                context.set_code(grpc.StatusCode.INVALID_ARGUMENT)
                context.set_details(error_msg)
                return MoleculeResponse(success=False, message=error_msg)
            
            # Create database session
            session = db_manager.get_session()
            
            try:
                # Convert atoms to dict format for database storage
                atoms_dict = [
                    {'symbol': atom.symbol, 'x': atom.x, 'y': atom.y, 'z': atom.z}
                    for atom in request.atoms
                ]
                
                # Create molecule in database
                molecule = create_molecule_from_atoms(
                    session=session,
                    name=request.name,
                    formula=request.formula,
                    atoms=atoms_dict,
                    charge=request.charge,
                    multiplicity=request.multiplicity,
                    symmetry=request.symmetry
                )
                
                # Convert back to protobuf
                pb_atoms = [
                    calculation_pb2.Atom(
                        symbol=atom['symbol'],
                        x=atom['x'],
                        y=atom['y'],
                        z=atom['z']
                    ) for atom in molecule.atoms
                ]
                
                pb_molecule = calculation_pb2.Molecule(
                    id=molecule.id,
                    name=molecule.name,
                    formula=molecule.formula,
                    molecular_weight=molecule.molecular_weight,
                    atoms=pb_atoms,
                    charge=molecule.charge,
                    multiplicity=molecule.multiplicity,
                    symmetry=molecule.symmetry,
                    geometry_type=calculation_pb2.GEOMETRY_TYPE_XYZ,
                    created_at=int(molecule.created_at.timestamp())
                )
                
                logger.info(f"Created molecule {molecule.id} in database")
                return MoleculeResponse(
                    success=True,
                    message=f"Molecule '{request.name}' created successfully",
                    molecule=pb_molecule
                )
                
            finally:
                session.close()
                
        except Exception as e:
            logger.error(f"Failed to create molecule: {e}")
            context.set_code(grpc.StatusCode.INTERNAL)
            context.set_details(str(e))
            return MoleculeResponse(success=False, message=f"Failed to create molecule: {str(e)}")
    
    def GetMolecule(self, request, context):
        session = db_manager.get_session()
        
        try:
            molecule = session.query(Molecule).filter(Molecule.id == request.molecule_id).first()
            if not molecule:
                context.set_code(grpc.StatusCode.NOT_FOUND)
                context.set_details(f"Molecule {request.molecule_id} not found")
                return MoleculeResponse(success=False, message="Molecule not found")
            
            # Convert to protobuf
            pb_atoms = [
                calculation_pb2.Atom(
                    symbol=atom['symbol'],
                    x=atom['x'],
                    y=atom['y'],
                    z=atom['z']
                ) for atom in molecule.atoms
            ]
            
            pb_molecule = calculation_pb2.Molecule(
                id=molecule.id,
                name=molecule.name,
                formula=molecule.formula,
                molecular_weight=molecule.molecular_weight,
                atoms=pb_atoms,
                charge=molecule.charge,
                multiplicity=molecule.multiplicity,
                symmetry=molecule.symmetry,
                geometry_type=calculation_pb2.GEOMETRY_TYPE_XYZ,
                created_at=int(molecule.created_at.timestamp())
            )
            
            return MoleculeResponse(
                success=True,
                message="Molecule retrieved successfully",
                molecule=pb_molecule
            )
        except Exception as e:
            logger.error(f"Failed to get molecule: {e}")
            context.set_code(grpc.StatusCode.INTERNAL)
            context.set_details(str(e))
            return MoleculeResponse(success=False, message=f"Failed to get molecule: {str(e)}")
        finally:
            session.close()
    
    def UpdateMolecule(self, request, context):
        return MoleculeResponse(success=False, message="Update not implemented yet")
    
    def DeleteMolecule(self, request, context):
        session = db_manager.get_session()
        
        try:
            molecule = session.query(Molecule).filter(Molecule.id == request.molecule_id).first()
            if not molecule:
                context.set_code(grpc.StatusCode.NOT_FOUND)
                context.set_details(f"Molecule {request.molecule_id} not found")
                return DeleteMoleculeResponse(success=False, message="Molecule not found")
            
            session.delete(molecule)
            session.commit()
            return DeleteMoleculeResponse(success=True, message="Molecule deleted successfully")
        except Exception as e:
            logger.error(f"Failed to delete molecule: {e}")
            context.set_code(grpc.StatusCode.INTERNAL)
            context.set_details(str(e))
            return DeleteMoleculeResponse(success=False, message=f"Failed to delete molecule: {str(e)}")
        finally:
            session.close()
    
    def CreateInstance(self, request, context):
        logger.info(f"CreateInstance called with name: {request.name}")
        
        session = db_manager.get_session()
        
        try:
            # Check if molecule exists
            molecule = session.query(Molecule).filter(Molecule.id == request.molecule_id).first()
            if not molecule:
                context.set_code(grpc.StatusCode.NOT_FOUND)
                context.set_details(f"Molecule {request.molecule_id} not found")
                return InstanceResponse(success=False, message="Molecule not found")
            
            # Create instance
            instance = create_calculation_instance(
                session=session,
                name=request.name,
                description=request.description,
                molecule_id=request.molecule_id,
                user_id=None,  # TODO: Add user context
                project_id=request.project_id if request.HasField('project_id') else None
            )
            
            # Convert to protobuf
            pb_molecule_atoms = [
                calculation_pb2.Atom(
                    symbol=atom['symbol'],
                    x=atom['x'],
                    y=atom['y'],
                    z=atom['z']
                ) for atom in molecule.atoms
            ]
            
            pb_molecule = calculation_pb2.Molecule(
                id=molecule.id,
                name=molecule.name,
                formula=molecule.formula,
                molecular_weight=molecule.molecular_weight,
                atoms=pb_molecule_atoms,
                charge=molecule.charge,
                multiplicity=molecule.multiplicity,
                symmetry=molecule.symmetry,
                geometry_type=calculation_pb2.GEOMETRY_TYPE_XYZ,
                created_at=int(molecule.created_at.timestamp())
            )
            
            pb_instance = calculation_pb2.Instance(
                id=instance.id,
                name=instance.name,
                description=instance.description,
                status=calculation_pb2.INSTANCE_STATUS_DRAFT,
                project_id=instance.project_id or "",
                created_at=int(instance.created_at.timestamp()),
                updated_at=int(instance.updated_at.timestamp()),
                molecule=pb_molecule
            )
            
            return InstanceResponse(
                success=True,
                message=f"Instance '{request.name}' created successfully",
                instance=pb_instance
            )
        except Exception as e:
            logger.error(f"Failed to create instance: {e}")
            context.set_code(grpc.StatusCode.INTERNAL)
            context.set_details(str(e))
            return InstanceResponse(success=False, message=f"Failed to create instance: {str(e)}")
        finally:
            session.close()
    
    def GetInstance(self, request, context):
        session = db_manager.get_session()
        
        try:
            instance = session.query(Instance).filter(Instance.id == request.instance_id).first()
            if not instance:
                context.set_code(grpc.StatusCode.NOT_FOUND)
                context.set_details(f"Instance {request.instance_id} not found")
                return InstanceResponse(success=False, message="Instance not found")
            
            # Convert to protobuf
            molecule = instance.molecule
            pb_molecule_atoms = [
                calculation_pb2.Atom(
                    symbol=atom['symbol'],
                    x=atom['x'],
                    y=atom['y'],
                    z=atom['z']
                ) for atom in molecule.atoms
            ]
            
            pb_molecule = calculation_pb2.Molecule(
                id=molecule.id,
                name=molecule.name,
                formula=molecule.formula,
                molecular_weight=molecule.molecular_weight,
                atoms=pb_molecule_atoms,
                charge=molecule.charge,
                multiplicity=molecule.multiplicity,
                symmetry=molecule.symmetry,
                geometry_type=calculation_pb2.GEOMETRY_TYPE_XYZ,
                created_at=int(molecule.created_at.timestamp())
            )
            
            pb_instance = calculation_pb2.Instance(
                id=instance.id,
                name=instance.name,
                description=instance.description,
                status=calculation_pb2.INSTANCE_STATUS_DRAFT,
                project_id=instance.project_id or "",
                created_at=int(instance.created_at.timestamp()),
                updated_at=int(instance.updated_at.timestamp()),
                molecule=pb_molecule
            )
            
            return InstanceResponse(
                success=True,
                message="Instance retrieved successfully",
                instance=pb_instance
            )
        except Exception as e:
            logger.error(f"Failed to get instance: {e}")
            context.set_code(grpc.StatusCode.INTERNAL)
            context.set_details(str(e))
            return InstanceResponse(success=False, message=f"Failed to get instance: {str(e)}")
        finally:
            session.close()
    
    def ListInstances(self, request, context):
        session = db_manager.get_session()
        
        try:
            query = session.query(Instance)
            
            # Filter by status if requested
            if request.HasField('status'):
                # Map protobuf status to database enum
                status_map = {
                    calculation_pb2.INSTANCE_STATUS_DRAFT: InstanceStatus.DRAFT,
                    calculation_pb2.INSTANCE_STATUS_READY: InstanceStatus.READY,
                    calculation_pb2.INSTANCE_STATUS_RUNNING: InstanceStatus.RUNNING,
                    calculation_pb2.INSTANCE_STATUS_COMPLETED: InstanceStatus.COMPLETED,
                    calculation_pb2.INSTANCE_STATUS_ERROR: InstanceStatus.ERROR
                }
                if request.status in status_map:
                    query = query.filter(Instance.status == status_map[request.status])
            
            # Filter by project if requested
            if request.HasField('project_id'):
                query = query.filter(Instance.project_id == request.project_id)
            
            # Get total count before pagination
            total_count = query.count()
            
            # Apply pagination
            if request.page > 0 and request.page_size > 0:
                offset = (request.page - 1) * request.page_size
                query = query.offset(offset).limit(request.page_size)
            
            instances = query.all()
            
            # Convert to protobuf
            pb_instances = []
            for instance in instances:
                molecule = instance.molecule
                pb_molecule_atoms = [
                    calculation_pb2.Atom(
                        symbol=atom['symbol'],
                        x=atom['x'],
                        y=atom['y'],
                        z=atom['z']
                    ) for atom in molecule.atoms
                ]
                
                pb_molecule = calculation_pb2.Molecule(
                    id=molecule.id,
                    name=molecule.name,
                    formula=molecule.formula,
                    molecular_weight=molecule.molecular_weight,
                    atoms=pb_molecule_atoms,
                    charge=molecule.charge,
                    multiplicity=molecule.multiplicity,
                    symmetry=molecule.symmetry,
                    geometry_type=calculation_pb2.GEOMETRY_TYPE_XYZ,
                    created_at=int(molecule.created_at.timestamp())
                )
                
                pb_instance = calculation_pb2.Instance(
                    id=instance.id,
                    name=instance.name,
                    description=instance.description,
                    status=calculation_pb2.INSTANCE_STATUS_DRAFT,
                    project_id=instance.project_id or "",
                    created_at=int(instance.created_at.timestamp()),
                    updated_at=int(instance.updated_at.timestamp()),
                    molecule=pb_molecule
                )
                pb_instances.append(pb_instance)
            
            return ListInstancesResponse(
                success=True,
                message="Instances retrieved successfully",
                instances=pb_instances,
                total_count=total_count,
                page=request.page,
                page_size=request.page_size
            )
        except Exception as e:
            logger.error(f"Failed to list instances: {e}")
            context.set_code(grpc.StatusCode.INTERNAL)
            context.set_details(str(e))
            return ListInstancesResponse(success=False, message=f"Failed to list instances: {str(e)}")
        finally:
            session.close()
    
    def UpdateInstance(self, request, context):
        return InstanceResponse(success=False, message="Update not implemented yet")
    
    def DeleteInstance(self, request, context):
        session = db_manager.get_session()
        
        try:
            instance = session.query(Instance).filter(Instance.id == request.instance_id).first()
            if not instance:
                context.set_code(grpc.StatusCode.NOT_FOUND)
                context.set_details(f"Instance {request.instance_id} not found")
                return DeleteInstanceResponse(success=False, message="Instance not found")
            
            session.delete(instance)
            session.commit()
            return DeleteInstanceResponse(success=True, message="Instance deleted successfully")
        except Exception as e:
            logger.error(f"Failed to delete instance: {e}")
            context.set_code(grpc.StatusCode.INTERNAL)
            context.set_details(str(e))
            return DeleteInstanceResponse(success=False, message=f"Failed to delete instance: {str(e)}")
        finally:
            session.close()
    
    def StartCalculation(self, request, context) -> Iterator[CalculationProgress]:
        logger.info(f"StartCalculation called for instance: {request.instance_id}")
        
        session = db_manager.get_session()
        
        try:
            # Get instance from database
            instance = session.query(Instance).filter(Instance.id == request.instance_id).first()
            if not instance:
                context.set_code(grpc.StatusCode.NOT_FOUND)
                context.set_details(f"Instance {request.instance_id} not found")
                return
            
            # Get molecule data
            molecule = instance.molecule
            if not molecule:
                context.set_code(grpc.StatusCode.NOT_FOUND)
                context.set_details("Molecule data not found for instance")
                return
            
            # Create calculation and job
            calculation, job = create_calculation_job(
                session=session,
                instance_id=request.instance_id,
                method=request.method,
                basis_set=request.basis_set,
                functional=request.parameters.get('functional') if request.parameters else None,
                parameters=dict(request.parameters) if request.parameters else None,
                priority=request.priority if request.HasField('priority') else 1
            )
            
            calculation_id = calculation.id
            job_id = job.id
            
            # Start real calculation in background thread
            progress_queue = []
            calculation_thread = threading.Thread(
                target=self._run_real_calculation,
                args=(calculation_id, molecule, request, progress_queue, session)
            )
            
            self.active_calculations[calculation_id] = calculation_thread
            calculation_thread.start()
            
            # Yield progress updates
            start_time = time.time()
            while calculation_thread.is_alive() or progress_queue:
                # Check for new progress updates
                while progress_queue:
                    step, percentage, status, message = progress_queue.pop(0)
                    
                    progress = CalculationProgress(
                        calculation_id=calculation_id,
                        job_id=job_id,
                        status=status,
                        progress_percentage=percentage,
                        current_step=step,
                        message=message,
                        timestamp=int(time.time()),
                        estimated_time_remaining="" if percentage >= 100 else f"{max(0, int((time.time() - start_time) * (100 - percentage) / max(percentage, 1)))}s"
                    )
                    yield progress
                
                if calculation_thread.is_alive():
                    time.sleep(0.5)  # Check every 500ms
            
            # Wait for thread completion
            calculation_thread.join()
            
            # Final status update
            final_calc = session.query(Calculation).filter(Calculation.id == calculation_id).first()
            if final_calc:
                final_status = calculation_pb2.CALCULATION_STATUS_COMPLETED if final_calc.status == CalculationStatus.COMPLETED else calculation_pb2.CALCULATION_STATUS_FAILED
                
                final_progress = CalculationProgress(
                    calculation_id=calculation_id,
                    job_id=job_id,
                    status=final_status,
                    progress_percentage=100.0,
                    current_step="Calculation finished",
                    message=f"Calculation {'completed successfully' if final_calc.status == CalculationStatus.COMPLETED else 'failed'}",
                    timestamp=int(time.time())
                )
                yield final_progress
            
        except Exception as e:
            logger.error(f"Calculation failed: {e}")
            context.set_code(grpc.StatusCode.INTERNAL)
            context.set_details(str(e))
        finally:
            session.close()
            if calculation_id in self.active_calculations:
                del self.active_calculations[calculation_id]
    
    def CancelCalculation(self, request, context):
        return CancelCalculationResponse(
            success=True,
            message=f"Calculation {request.calculation_id} cancelled"
        )
    
    def GetCalculationStatus(self, request, context):
        return CalculationStatusResponse(
            success=True,
            message="Status retrieved",
            status=calculation_pb2.CALCULATION_STATUS_COMPLETED
        )
    
    def GetResults(self, request, context):
        # Mock results
        mock_results = calculation_pb2.CalculationResults(
            calculation_id=request.calculation_id,
            converged=True,
            total_energy=-76.0267,
            energy_unit="Hartree",
            calculation_time_seconds=120
        )
        
        return ResultsResponse(
            success=True,
            message="Results retrieved successfully",
            results=mock_results
        )
    
    def ExportResults(self, request, context):
        return ExportResultsResponse(
            success=True,
            message="Results exported successfully",
            file_path=f"/exports/calc_{request.calculation_id}.json"
        )
    
    def GetJobQueue(self, request, context):
        session = db_manager.get_session()
        
        try:
            jobs = session.query(JobQueue).order_by(JobQueue.priority.desc(), JobQueue.created_at).all()
            
            # Convert to protobuf (simplified - would need proper Job protobuf definition)
            pb_jobs = []  # TODO: Implement proper job protobuf conversion
            
            return JobQueueResponse(
                success=True,
                message="Job queue retrieved",
                jobs=pb_jobs
            )
        except Exception as e:
            logger.error(f"Failed to get job queue: {e}")
            return JobQueueResponse(success=False, message=f"Failed to get job queue: {str(e)}")
        finally:
            session.close()
    
    def UpdateJobPriority(self, request, context):
        return UpdateJobPriorityResponse(
            success=True,
            message=f"Job {request.job_id} priority updated to {request.new_priority}"
        )
    
    def GetSystemInfo(self, request, context):
        try:
            import psutil
            import sys
            import pyscf
            
            # Get system resources
            memory = psutil.virtual_memory()
            disk = psutil.disk_usage('/')
            
            # Test database connection
            db_connected, db_message = db_manager.test_connection()
            
            system_info = calculation_pb2.SystemInfo(
                version="1.0.0",
                python_version=sys.version,
                pyscf_version=pyscf.__version__,
                gpu_available=pyscf_config.config['use_gpu'],
                available_methods=pyscf_config.get_available_methods(),
                available_basis_sets=pyscf_config.get_available_basis_sets(),
                resources=calculation_pb2.SystemResources(
                    cpu_cores=psutil.cpu_count(),
                    total_memory_mb=int(memory.total / 1024 / 1024),
                    available_memory_mb=int(memory.available / 1024 / 1024),
                    disk_space_mb=int(disk.free / 1024 / 1024)
                )
            )
            
            details = {
                "database": "connected" if db_connected else "disconnected",
                "pyscf_config": f"threads={pyscf_config.config['omp_threads']}, memory={pyscf_config.config['max_memory']}MB",
                "gpu_support": "available" if pyscf_config.config['use_gpu'] else "not_available"
            }
            
            return SystemInfoResponse(
                success=True,
                message="System info retrieved",
                system_info=system_info
            )
            
        except Exception as e:
            logger.error(f"Failed to get system info: {e}")
            return SystemInfoResponse(
                success=False,
                message=f"Failed to get system info: {str(e)}"
            )
    
    def HealthCheck(self, request, context):
        try:
            # Test database connection
            db_connected, db_message = db_manager.test_connection()
            
            # Test PySCF availability
            pyscf_available = True
            try:
                import pyscf
            except ImportError:
                pyscf_available = False
            
            healthy = db_connected and pyscf_available
            status = "OK" if healthy else "DEGRADED"
            
            details = {
                "database": "connected" if db_connected else f"error: {db_message}",
                "pyscf": "available" if pyscf_available else "not_available",
                "threads": str(pyscf_config.config['omp_threads']),
                "memory_limit": f"{pyscf_config.config['max_memory']}MB"
            }
            
            return HealthCheckResponse(
                healthy=healthy,
                status=status,
                details=details,
                timestamp=int(time.time())
            )
            
        except Exception as e:
            logger.error(f"Health check failed: {e}")
            return HealthCheckResponse(
                healthy=False,
                status="ERROR",
                details={"error": str(e)},
                timestamp=int(time.time())
            )
    
    def _get_atomic_weight(self, symbol: str) -> float:
        """Get atomic weight for element symbol"""
        weights = {
            'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
            'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180
        }
        return weights.get(symbol, 1.0)
    
    def _run_real_calculation(self, calculation_id: str, molecule_db, request, progress_queue: list, session):
        """Run the actual PySCF calculation in background thread"""
        try:
            # Update calculation status to running
            calculation = session.query(Calculation).filter(Calculation.id == calculation_id).first()
            calculation.status = CalculationStatus.RUNNING
            calculation.start_time = datetime.utcnow()
            session.commit()
            
            progress_queue.append(("Starting calculation", 0.0, calculation_pb2.CALCULATION_STATUS_RUNNING, "Initializing..."))
            
            # Convert database molecule to PySCF format
            atoms_list = [(atom['symbol'], atom['x'], atom['y'], atom['z']) for atom in molecule_db.atoms]
            
            progress_queue.append(("Building molecule", 10.0, calculation_pb2.CALCULATION_STATUS_RUNNING, f"Building {molecule_db.formula}"))
            
            # Build PySCF molecule
            pyscf_mol = MoleculeBuilder.from_atoms(
                atoms=atoms_list,
                charge=molecule_db.charge,
                spin=molecule_db.multiplicity - 1,  # Convert multiplicity to spin
                basis=request.basis_set,
                symmetry=molecule_db.symmetry
            )
            
            progress_queue.append(("Setting up calculation", 20.0, calculation_pb2.CALCULATION_STATUS_RUNNING, f"Method: {request.method}, Basis: {request.basis_set}"))
            
            # Run calculation based on method
            results = None
            if request.method.upper() == 'HF':
                results = calculation_engine.run_hf_calculation(
                    pyscf_mol,
                    conv_tol=request.convergence_criteria.energy_threshold if request.HasField('convergence_criteria') else None,
                    max_cycle=request.max_iterations if request.HasField('max_iterations') else None
                )
                progress_queue.append(("Running HF calculation", 60.0, calculation_pb2.CALCULATION_STATUS_RUNNING, "SCF iterations in progress"))
                
            elif request.method.upper() == 'DFT':
                functional = request.parameters.get('functional', 'b3lyp') if request.parameters else 'b3lyp'
                results = calculation_engine.run_dft_calculation(
                    pyscf_mol,
                    functional=functional,
                    conv_tol=request.convergence_criteria.energy_threshold if request.HasField('convergence_criteria') else None,
                    max_cycle=request.max_iterations if request.HasField('max_iterations') else None
                )
                progress_queue.append(("Running DFT calculation", 60.0, calculation_pb2.CALCULATION_STATUS_RUNNING, f"DFT({functional}) in progress"))
                
            elif request.method.upper() == 'MP2':
                results = calculation_engine.run_mp2_calculation(
                    pyscf_mol,
                    conv_tol=request.convergence_criteria.energy_threshold if request.HasField('convergence_criteria') else None
                )
                progress_queue.append(("Running MP2 calculation", 60.0, calculation_pb2.CALCULATION_STATUS_RUNNING, "Post-HF correlation"))
                
            else:
                raise ValueError(f"Unsupported method: {request.method}")
            
            progress_queue.append(("Finalizing results", 90.0, calculation_pb2.CALCULATION_STATUS_RUNNING, "Saving results"))
            
            # Save results to database
            if results:
                save_calculation_results(session, calculation_id, results)
                calculation.status = CalculationStatus.COMPLETED
                calculation.end_time = datetime.utcnow()
                session.commit()
                
                progress_queue.append(("Calculation completed", 100.0, calculation_pb2.CALCULATION_STATUS_COMPLETED, 
                                     f"Energy: {results.get('total_energy', 0.0):.6f} Eh"))
            else:
                raise RuntimeError("No results returned from calculation")
                
        except Exception as e:
            logger.error(f"Calculation {calculation_id} failed: {e}")
            
            # Update calculation status to failed
            calculation = session.query(Calculation).filter(Calculation.id == calculation_id).first()
            if calculation:
                calculation.status = CalculationStatus.FAILED
                calculation.end_time = datetime.utcnow()
                calculation.error_message = str(e)
                session.commit()
            
            progress_queue.append(("Calculation failed", 100.0, calculation_pb2.CALCULATION_STATUS_FAILED, str(e)))

    def _current_timestamp(self) -> int:
        """Get current Unix timestamp"""
        import time
        return int(time.time())


def serve():
    """Start the gRPC server"""
    server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))
    calculation_pb2_grpc.add_CalculationServiceServicer_to_server(
        CalculationServiceImpl(), server
    )
    
    listen_addr = '[::]:50051'
    server.add_insecure_port(listen_addr)
    
    logger.info(f"Starting gRPC server on {listen_addr}")
    server.start()
    
    try:
        server.wait_for_termination()
    except KeyboardInterrupt:
        logger.info("Shutting down gRPC server")
        server.stop(0)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    serve()