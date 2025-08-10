from pydantic import BaseModel, Field, validator
from typing import Optional, Literal
from enum import Enum


class SearchTypeEnum(str, Enum):
    name = "name"
    cid = "cid"
    formula = "formula"


class SolventMethodEnum(str, Enum):
    none = "none"
    pcm = "pcm"


class PubChemSearchModel(BaseModel):
    query: str = Field(..., min_length=1, description="Search query for PubChem")
    search_type: SearchTypeEnum = Field(default=SearchTypeEnum.name, description="Type of search: name, cid, or formula")


class SMILESConvertModel(BaseModel):
    smiles: str = Field(..., min_length=1, description="SMILES string to convert to XYZ format")


class PubChemValidateModel(BaseModel):
    xyz: str = Field(..., min_length=1, description="XYZ string to validate")


class QuantumCalculateModel(BaseModel):
    xyz: str = Field(..., min_length=1, description="XYZ molecular structure data")
    calculation_method: Literal["DFT"] = Field(default="DFT", description="Calculation method (only DFT supported)")
    basis_function: str = Field(default="6-31G(d)", description="Basis function for calculation")
    exchange_correlation: str = Field(default="B3LYP", description="Exchange-correlation functional")
    charges: int = Field(default=0, ge=-10, le=10, description="Molecular charge (-10 to +10)")
    spin_multiplicity: int = Field(default=1, ge=1, le=10, description="Spin multiplicity (1-10)")
    solvent_method: SolventMethodEnum = Field(default=SolventMethodEnum.none, description="Solvent method")
    solvent: str = Field(default="-", description="Solvent type")
    molecule_name: str = Field(default="Unnamed Calculation", min_length=1, max_length=100, description="Name for the calculation")
    cpu_cores: Optional[int] = Field(default=None, ge=1, le=32, description="Number of CPU cores (1-32)")
    memory_mb: Optional[int] = Field(default=None, ge=512, le=32768, description="Memory in MB (512-32768)")

    @validator('xyz')
    def validate_xyz(cls, v):
        lines = v.strip().split('\n')
        if len(lines) < 3:
            raise ValueError('XYZ data must contain at least 3 lines (atom count, comment, and one atom)')
        try:
            atom_count = int(lines[0].strip())
            if atom_count < 1:
                raise ValueError('Atom count must be at least 1')
            if len(lines) < atom_count + 2:
                raise ValueError(f'XYZ data should have {atom_count + 2} lines for {atom_count} atoms')
        except ValueError as e:
            if "invalid literal for int()" in str(e):
                raise ValueError('First line of XYZ data must be a valid integer representing atom count')
            raise e
        return v.strip()

    @validator('molecule_name')
    def validate_molecule_name(cls, v):
        if not v or not v.strip():
            return "Unnamed Calculation"
        return v.strip()


class CalculationUpdateModel(BaseModel):
    name: str = Field(..., min_length=1, max_length=100, description="Updated name for the calculation")

    @validator('name')
    def validate_name(cls, v):
        if not v or not v.strip():
            raise ValueError('Calculation name cannot be empty')
        return v.strip()