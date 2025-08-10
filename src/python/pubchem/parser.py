"""XYZ format parser and converter for molecular structures."""

from typing import List, Dict, Any
import logging

logger = logging.getLogger(__name__)


def atoms_to_xyz(atoms: List[Dict[str, Any]], title: str = "Molecule") -> str:
    """Convert list of atoms to XYZ format string.
    
    Args:
        atoms: List of atom dictionaries with 'element', 'x', 'y', 'z' keys.
        title: Title/comment line for the XYZ file.
        
    Returns:
        XYZ format string.
        
    Raises:
        ValueError: If atom data is invalid.
    """
    if not atoms:
        raise ValueError("No atoms provided to convert to XYZ format.")

    # Validate atom data
    for i, atom in enumerate(atoms):
        if not all(k in atom for k in ['element', 'x', 'y', 'z']):
            raise ValueError(f"Atom {i} is missing required coordinate keys.")
        if not isinstance(atom.get('element'), str) or not all(isinstance(atom.get(c), (int, float)) for c in ['x', 'y', 'z']):
            raise ValueError(f"Atom {i} has invalid data types.")

    # Build XYZ string
    lines = [str(len(atoms)), title]
    for atom in atoms:
        lines.append(f"{atom['element']:<3} {atom['x']:>10.4f} {atom['y']:>10.4f} {atom['z']:>10.4f}")

    return '\n'.join(lines)


def validate_xyz(xyz_string: str) -> Dict[str, Any]:
    """Validate XYZ format string.
    
    Returns:
        Dictionary with validation results:
        - 'valid': bool
        - 'error': str (if invalid)
        - 'num_atoms': int (if valid)
        - 'atoms': List[Dict] (if valid)
    """
    try:
        lines = xyz_string.strip().split('\n')
        
        if len(lines) < 2:
            return {'valid': False, 'error': 'XYZ string must have at least 2 lines (count and title).'}
        
        num_atoms = int(lines[0].strip())
        if num_atoms <= 0:
            return {'valid': False, 'error': 'Number of atoms must be a positive integer.'}

        expected_lines = num_atoms + 2
        if len(lines) != expected_lines:
            return {'valid': False, 'error': f'Expected {expected_lines} lines for {num_atoms} atoms, but got {len(lines)}.'}
        
        atoms = []
        for i, line_str in enumerate(lines[2:]):
            parts = line_str.strip().split()
            if len(parts) != 4:
                return {'valid': False, 'error': f'Atom line {i+3} should have 4 parts (element, x, y, z).'}
            
            element, x, y, z = parts
            atoms.append({'element': element, 'x': float(x), 'y': float(y), 'z': float(z)})
        
        return {'valid': True, 'num_atoms': num_atoms, 'atoms': atoms, 'title': lines[1].strip()}

    except (ValueError, IndexError) as e:
        return {'valid': False, 'error': f'Parsing failed. Ensure correct numeric formats. Error: {e}'}
    except Exception as e:
        return {'valid': False, 'error': f'An unexpected validation error occurred: {e}'}


def format_compound_title(compound_data: Any, query: str) -> str:
    """Format a title for an XYZ file from a CompoundData object.
    
    Args:
        compound_data: A CompoundData object (or similar structure).
        query: The original user query string.
    
    Returns:
        Formatted title string for the XYZ file.
    """
    cid = getattr(compound_data, 'cid', 'N/A')
    formula = getattr(compound_data, 'molecular_formula', '')
    name = getattr(compound_data, 'iupac_name', query) or query

    # Truncate long names
    if len(name) > 60:
        name = name[:57] + "..."
    
    if formula and formula != 'Unknown':
        return f"{name} ({formula}) - CID: {cid}"
    else:
        return f"{name} - CID: {cid}"