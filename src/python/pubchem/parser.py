"""XYZ format parser and converter for molecular structures."""

from typing import List, Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)


class XYZParser:
    """Parser for converting molecular structure data to XYZ format."""
    
    @staticmethod
    def atoms_to_xyz(atoms: List[Dict[str, Any]], title: str = "Molecule") -> str:
        """Convert list of atoms to XYZ format string.
        
        Args:
            atoms: List of atom dictionaries with 'element', 'x', 'y', 'z' keys
            title: Title/comment line for the XYZ file
            
        Returns:
            XYZ format string
            
        Raises:
            ValueError: If atom data is invalid
        """
        if not atoms:
            raise ValueError("No atoms provided")
        
        # Validate atom data
        for i, atom in enumerate(atoms):
            required_keys = ['element', 'x', 'y', 'z']
            missing_keys = [key for key in required_keys if key not in atom]
            if missing_keys:
                raise ValueError(f"Atom {i} missing keys: {missing_keys}")
            
            # Check coordinate types
            for coord in ['x', 'y', 'z']:
                if not isinstance(atom[coord], (int, float)):
                    raise ValueError(f"Atom {i} has non-numeric {coord}: {atom[coord]}")
        
        # Build XYZ string
        lines = []
        lines.append(str(len(atoms)))  # Number of atoms
        lines.append(title)  # Comment line
        
        # Atom lines with compact formatting
        for atom in atoms:
            element = atom['element']
            x, y, z = atom['x'], atom['y'], atom['z']
            lines.append(f"{element:<2s} {x:8.4f} {y:8.4f} {z:8.4f}")

        return '\n'.join(lines)
    
    @staticmethod
    def validate_xyz(xyz_string: str) -> Dict[str, Any]:
        """Validate XYZ format string.
        
        Args:
            xyz_string: XYZ format string
            
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
                return {'valid': False, 'error': 'XYZ string too short'}
            
            # Parse number of atoms
            try:
                num_atoms = int(lines[0].strip())
            except ValueError:
                return {'valid': False, 'error': 'First line must be number of atoms'}
            
            if num_atoms <= 0:
                return {'valid': False, 'error': 'Number of atoms must be positive'}
            
            # Check if we have enough lines
            expected_lines = num_atoms + 2  # atoms + count line + comment line
            if len(lines) < expected_lines:
                return {'valid': False, 'error': f'Expected {expected_lines} lines, got {len(lines)}'}
            
            # Parse atoms
            atoms = []
            for i in range(2, 2 + num_atoms):  # Skip count and comment lines
                line = lines[i].strip()
                if not line:
                    return {'valid': False, 'error': f'Empty atom line at {i+1}'}
                
                parts = line.split()
                if len(parts) < 4:
                    return {'valid': False, 'error': f'Atom line {i+1} has insufficient data'}
                
                try:
                    element = parts[0]
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    
                    atoms.append({
                        'element': element,
                        'x': x,
                        'y': y,
                        'z': z
                    })
                except ValueError as e:
                    return {'valid': False, 'error': f'Invalid coordinates at line {i+1}: {e}'}
            
            return {
                'valid': True,
                'num_atoms': num_atoms,
                'atoms': atoms,
                'title': lines[1].strip() if len(lines) > 1 else ''
            }
            
        except Exception as e:
            return {'valid': False, 'error': f'Unexpected error: {e}'}
    
    @staticmethod
    def format_compound_title(compound_info: Dict[str, Any], query: str) -> str:
        """Format a title for XYZ file from compound information.
        
        Args:
            compound_info: Compound information dictionary
            query: Original query string
            
        Returns:
            Formatted title string
        """
        cid = compound_info.get('cid', 'Unknown')
        formula = compound_info.get('molecular_formula', '')
        iupac_name = compound_info.get('iupac_name', query)
        
        # Truncate long names
        if len(iupac_name) > 50:
            iupac_name = iupac_name[:47] + "..."
        
        if formula and formula != 'Unknown':
            return f"{iupac_name} ({formula}) - CID:{cid}"
        else:
            return f"{iupac_name} - CID:{cid}"