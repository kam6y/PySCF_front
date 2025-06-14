"""
SQLAlchemy models for PySCF database
"""
import uuid
from datetime import datetime
from typing import Optional, List
from sqlalchemy import Column, Integer, String, Text, Float, DateTime, Boolean, ForeignKey, JSON, Enum
from sqlalchemy.orm import relationship
from sqlalchemy.dialects.mysql import CHAR
import enum

from .connection import Base

# Enums
class InstanceStatus(enum.Enum):
    DRAFT = "draft"
    READY = "ready"
    RUNNING = "running"
    COMPLETED = "completed"
    ERROR = "error"

class CalculationStatus(enum.Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"

class JobStatus(enum.Enum):
    WAITING = "waiting"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"

# Database Models
class Molecule(Base):
    __tablename__ = "molecules"
    
    id = Column(CHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String(255), nullable=False)
    formula = Column(String(100), nullable=False)
    molecular_weight = Column(Float, nullable=True)
    charge = Column(Integer, default=0)
    multiplicity = Column(Integer, default=1)
    symmetry = Column(String(10), default='c1')
    geometry_type = Column(String(20), default='xyz')
    atoms = Column(JSON, nullable=False)  # List of atom dictionaries
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    instances = relationship("Instance", back_populates="molecule", cascade="all, delete-orphan")

class Instance(Base):
    __tablename__ = "instances"
    
    id = Column(CHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    status = Column(Enum(InstanceStatus), default=InstanceStatus.DRAFT)
    user_id = Column(String(100), nullable=True)
    project_id = Column(String(100), nullable=True)
    molecule_id = Column(CHAR(36), ForeignKey('molecules.id'), nullable=False)
    settings = Column(JSON, nullable=True)  # Additional settings
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    molecule = relationship("Molecule", back_populates="instances")
    calculations = relationship("Calculation", back_populates="instance", cascade="all, delete-orphan")

class Calculation(Base):
    __tablename__ = "calculations"
    
    id = Column(CHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    instance_id = Column(CHAR(36), ForeignKey('instances.id'), nullable=False)
    method = Column(String(50), nullable=False)
    basis_set = Column(String(50), nullable=False)
    functional = Column(String(50), nullable=True)  # For DFT calculations
    parameters = Column(JSON, nullable=True)
    convergence_criteria = Column(JSON, nullable=True)
    max_iterations = Column(Integer, default=50)
    start_time = Column(DateTime, nullable=True)
    end_time = Column(DateTime, nullable=True)
    status = Column(Enum(CalculationStatus), default=CalculationStatus.PENDING)
    error_message = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    
    # Relationships
    instance = relationship("Instance", back_populates="calculations")
    results = relationship("Result", back_populates="calculation", cascade="all, delete-orphan")
    job = relationship("JobQueue", back_populates="calculation", uselist=False)

class Result(Base):
    __tablename__ = "results"
    
    id = Column(CHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    calculation_id = Column(CHAR(36), ForeignKey('calculations.id'), nullable=False)
    converged = Column(Boolean, default=False)
    total_energy = Column(Float, nullable=True)
    energy_unit = Column(String(20), default='Hartree')
    calculation_time_seconds = Column(Integer, nullable=True)
    iterations = Column(Integer, nullable=True)
    
    # Energy components
    nuclear_repulsion_energy = Column(Float, nullable=True)
    electronic_energy = Column(Float, nullable=True)
    correlation_energy = Column(Float, nullable=True)
    
    # Orbital information
    homo_energy = Column(Float, nullable=True)
    lumo_energy = Column(Float, nullable=True)
    homo_lumo_gap = Column(Float, nullable=True)
    orbital_energies = Column(JSON, nullable=True)
    occupations = Column(JSON, nullable=True)
    
    # Molecular properties
    dipole_moment = Column(JSON, nullable=True)
    quadrupole_moment = Column(JSON, nullable=True)
    atomic_charges = Column(JSON, nullable=True)
    
    # Frequencies (if applicable)
    frequencies = Column(JSON, nullable=True)
    intensities = Column(JSON, nullable=True)
    
    # Raw data (binary)
    molecular_orbitals = Column(Text, nullable=True)  # Base64 encoded
    normal_modes = Column(Text, nullable=True)  # Base64 encoded
    
    # Additional data
    additional_data = Column(JSON, nullable=True)
    
    created_at = Column(DateTime, default=datetime.utcnow)
    
    # Relationships
    calculation = relationship("Calculation", back_populates="results")

class JobQueue(Base):
    __tablename__ = "job_queue"
    
    id = Column(CHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    calculation_id = Column(CHAR(36), ForeignKey('calculations.id'), nullable=False)
    priority = Column(Integer, default=1)
    status = Column(Enum(JobStatus), default=JobStatus.WAITING)
    assigned_worker = Column(String(100), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)
    error_message = Column(Text, nullable=True)
    progress_data = Column(JSON, nullable=True)
    
    # Relationships
    calculation = relationship("Calculation", back_populates="job")

# Helper functions for model operations
def create_molecule_from_atoms(session, name: str, formula: str, atoms: List[dict], 
                              charge: int = 0, multiplicity: int = 1, 
                              symmetry: str = 'c1') -> Molecule:
    """Create a molecule from atom list"""
    
    # Calculate molecular weight (simplified)
    atomic_weights = {
        'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
        'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
        'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
        'S': 32.066, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078
    }
    
    molecular_weight = sum(atomic_weights.get(atom['symbol'], 0.0) for atom in atoms)
    
    molecule = Molecule(
        name=name,
        formula=formula,
        molecular_weight=molecular_weight,
        charge=charge,
        multiplicity=multiplicity,
        symmetry=symmetry,
        atoms=atoms
    )
    
    session.add(molecule)
    session.commit()
    session.refresh(molecule)
    
    return molecule

def create_calculation_instance(session, name: str, description: str, 
                               molecule_id: str, user_id: Optional[str] = None,
                               project_id: Optional[str] = None) -> Instance:
    """Create a calculation instance"""
    
    instance = Instance(
        name=name,
        description=description,
        molecule_id=molecule_id,
        user_id=user_id,
        project_id=project_id,
        status=InstanceStatus.DRAFT
    )
    
    session.add(instance)
    session.commit()
    session.refresh(instance)
    
    return instance

def create_calculation_job(session, instance_id: str, method: str, basis_set: str,
                          functional: Optional[str] = None, parameters: Optional[dict] = None,
                          priority: int = 1) -> tuple[Calculation, JobQueue]:
    """Create a calculation and its job queue entry"""
    
    calculation = Calculation(
        instance_id=instance_id,
        method=method,
        basis_set=basis_set,
        functional=functional,
        parameters=parameters or {},
        status=CalculationStatus.PENDING
    )
    
    session.add(calculation)
    session.commit()
    session.refresh(calculation)
    
    job = JobQueue(
        calculation_id=calculation.id,
        priority=priority,
        status=JobStatus.WAITING
    )
    
    session.add(job)
    session.commit()
    session.refresh(job)
    
    return calculation, job

def save_calculation_results(session, calculation_id: str, results_data: dict) -> Result:
    """Save calculation results"""
    
    result = Result(
        calculation_id=calculation_id,
        converged=results_data.get('converged', False),
        total_energy=results_data.get('total_energy'),
        nuclear_repulsion_energy=results_data.get('nuclear_repulsion'),
        electronic_energy=results_data.get('electronic_energy'),
        correlation_energy=results_data.get('correlation_energy'),
        homo_energy=results_data.get('homo_energy'),
        lumo_energy=results_data.get('lumo_energy'),
        orbital_energies=results_data.get('orbital_energies'),
        occupations=results_data.get('occupations'),
        iterations=results_data.get('cycles', 0),
        additional_data=results_data.get('additional_data', {})
    )
    
    # Calculate HOMO-LUMO gap if both are available
    if result.homo_energy is not None and result.lumo_energy is not None:
        result.homo_lumo_gap = result.lumo_energy - result.homo_energy
    
    session.add(result)
    session.commit()
    session.refresh(result)
    
    return result