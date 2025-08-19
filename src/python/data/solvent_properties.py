"""Solvent properties data for quantum chemistry calculations."""

from typing import Dict

# Dielectric constants for common solvents
SOLVENT_DIELECTRIC: Dict[str, float] = {
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