"""Solvent effects setup for quantum chemistry calculations using PySCF."""

import logging
from typing import Union, Any
from pyscf import scf, dft

logger = logging.getLogger(__name__)


def setup_pcm_solvent(mf: Union[scf.hf.SCF, dft.rks.RKS], method: str, solvent: str) -> None:
    """Setup PCM solvent effects for a PySCF mean field object.
    
    Args:
        mf: PySCF mean field object (SCF or DFT)
        method: PCM method ('pcm', 'ief-pcm', 'c-pcm', 'cosmo', 'ssvpe')
        solvent: Solvent name or dielectric constant
    """
    # Set PCM method
    method_lower = method.lower()
    if method_lower == 'ief-pcm':
        mf.with_solvent.method = 'IEF-PCM'
    elif method_lower == 'c-pcm':
        mf.with_solvent.method = 'C-PCM'
    elif method_lower == 'cosmo':
        mf.with_solvent.method = 'COSMO'
    elif method_lower == 'ssvpe':
        mf.with_solvent.method = 'SS(V)PE'
    else:  # default to IEF-PCM for 'pcm'
        mf.with_solvent.method = 'IEF-PCM'
    
    # Set dielectric constant based on solvent
    solvent_dielectric = {
        'water': 78.3553,
        'dimethylsulfoxide': 46.826,
        'dmso': 46.826,  # alias
        'n,n-dimethylformamide': 37.219,
        'dmf': 37.219,  # alias
        'nitromethane': 36.562,
        'methanol': 32.613,
        'ethanol': 24.852,
        'acetone': 20.493,
        'dichloroethane': 10.125,
        'dichloromethane': 8.93,
        'tetrahydrofuran': 7.4297,
        'thf': 7.4297,  # alias
        'chlorobenzene': 5.6968,
        'chloroform': 4.7113,
        'diethylether': 4.2400,
        'toluene': 2.3741,
        'benzene': 2.2706,
        '1,4-dioxane': 2.2099,
        'dioxane': 2.2099,  # alias
        'cyclohexane': 2.0160
    }
                
    if solvent.lower() in solvent_dielectric:
        mf.with_solvent.eps = solvent_dielectric[solvent.lower()]
    else:
        # Try to parse as custom dielectric constant
        try:
            eps_value = float(solvent)
            if eps_value > 1.0:
                mf.with_solvent.eps = eps_value
            else:
                logger.warning(f"Invalid dielectric constant: {solvent}, using water (78.36)")
                mf.with_solvent.eps = 78.3553
        except ValueError:
            logger.warning(f"Unknown solvent: {solvent}, using water dielectric constant")
            mf.with_solvent.eps = 78.3553


def setup_ddcosmo_solvent(mf: Union[scf.hf.SCF, dft.rks.RKS], solvent: str) -> Union[scf.hf.SCF, dft.rks.RKS]:
    """Setup ddCOSMO solvent effects for a PySCF mean field object.
    
    Args:
        mf: PySCF mean field object (SCF or DFT)
        solvent: Solvent name or dielectric constant
        
    Returns:
        Modified mean field object with ddCOSMO solvent effects
    """
    # Apply ddCOSMO solvent model using PySCF's DDCOSMO method
    mf = mf.DDCOSMO()
    
    # Set dielectric constant (reuse the same dictionary as PCM)
    solvent_dielectric = {
        'water': 78.3553,
        'dimethylsulfoxide': 46.826,
        'dmso': 46.826,  # alias
        'n,n-dimethylformamide': 37.219,
        'dmf': 37.219,  # alias
        'nitromethane': 36.562,
        'methanol': 32.613,
        'ethanol': 24.852,
        'acetone': 20.493,
        'dichloroethane': 10.125,
        'dichloromethane': 8.93,
        'tetrahydrofuran': 7.4297,
        'thf': 7.4297,  # alias
        'chlorobenzene': 5.6968,
        'chloroform': 4.7113,
        'diethylether': 4.2400,
        'toluene': 2.3741,
        'benzene': 2.2706,
        '1,4-dioxane': 2.2099,
        'dioxane': 2.2099,  # alias
        'cyclohexane': 2.0160
    }
    
    if solvent.lower() in solvent_dielectric:
        eps = solvent_dielectric[solvent.lower()]
    else:
        # Try to parse as custom dielectric constant
        try:
            eps_value = float(solvent)
            if eps_value > 1.0:
                eps = eps_value
            else:
                logger.warning(f"Invalid dielectric constant: {solvent}, using water (78.36)")
                eps = 78.3553
        except ValueError:
            logger.warning(f"Unknown solvent: {solvent}, using water dielectric constant")
            eps = 78.3553
    
    # Set ddCOSMO parameters
    mf.with_solvent.eps = eps
    
    # Optional: Set additional ddCOSMO parameters
    # mf.with_solvent.lebedev_order = 29  # default spherical grid
    # mf.with_solvent.lmax = 10  # default maximum angular momentum
    
    return mf


def setup_solvent_effects(mf: Union[scf.hf.SCF, dft.rks.RKS], method: str, solvent: str) -> Union[scf.hf.SCF, dft.rks.RKS]:
    """Setup solvent effects for a PySCF mean field object.
    
    Args:
        mf: PySCF mean field object (SCF or DFT)
        method: Solvent method ('none', 'ief-pcm', 'c-pcm', 'cosmo', 'ssvpe', 'ddcosmo')
        solvent: Solvent name or dielectric constant
        
    Returns:
        Modified mean field object with solvent effects
    """
    if method == 'none' or solvent == '-':
        return mf
    elif method.lower() in ['ief-pcm', 'c-pcm', 'cosmo', 'ssvpe']:
        mf = mf.PCM()
        setup_pcm_solvent(mf, method, solvent)
        return mf
    elif method.lower() == 'ddcosmo':
        return setup_ddcosmo_solvent(mf, solvent)
    else:
        logger.warning(f"Unsupported solvent method: {method}, using no solvent")
        return mf