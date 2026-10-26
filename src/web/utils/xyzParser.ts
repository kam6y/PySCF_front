export interface Atom {
  element: string;
  x: number;
  y: number;
  z: number;
}

export interface ParsedXYZ {
  numAtoms: number;
  comment: string;
  atoms: Atom[];
}

export interface XYZValidationResult {
  isValid: boolean;
  error?: string;
  data?: ParsedXYZ;
}

/**
 * Parse XYZ coordinate data from string
 * XYZ format:
 * Line 1: Number of atoms
 * Line 2: Comment (optional)
 * Lines 3+: Element X Y Z
 */
export function parseXYZData(xyzString: string): XYZValidationResult {
  try {
    const lines = xyzString
      .trim()
      .split('\n')
      .filter(line => line.trim() !== '');

    if (lines.length < 2) {
      return {
        isValid: false,
        error:
          'XYZ data must have at least 2 lines (atom count and at least one atom)',
      };
    }

    // Parse number of atoms
    const numAtomsLine = lines[0].trim();
    const numAtoms = parseInt(numAtomsLine, 10);

    if (isNaN(numAtoms) || numAtoms <= 0) {
      return {
        isValid: false,
        error: 'First line must be a positive integer (number of atoms)',
      };
    }

    // Parse comment line (optional)
    let comment = '';
    let atomStartIndex = 1;

    // Check if second line looks like atom data or comment
    if (lines.length > 1) {
      const secondLineParts = lines[1].trim().split(/\s+/);

      // If second line has exactly 4 parts and parts 1-3 are numbers, it's atom data
      if (
        secondLineParts.length === 4 &&
        !isNaN(parseFloat(secondLineParts[1])) &&
        !isNaN(parseFloat(secondLineParts[2])) &&
        !isNaN(parseFloat(secondLineParts[3]))
      ) {
        comment = '';
        atomStartIndex = 1;
      } else {
        comment = lines[1].trim();
        atomStartIndex = 2;
      }
    }

    // Parse atoms
    const atoms: Atom[] = [];
    const atomLines = lines.slice(atomStartIndex);

    if (atomLines.length !== numAtoms) {
      return {
        isValid: false,
        error: `Expected ${numAtoms} atoms, but found ${atomLines.length} atom lines`,
      };
    }

    for (let i = 0; i < atomLines.length; i++) {
      const line = atomLines[i].trim();
      const parts = line.split(/\s+/);

      if (parts.length < 4) {
        return {
          isValid: false,
          error: `Line ${atomStartIndex + i + 1}: Expected format "Element X Y Z", got "${line}"`,
        };
      }

      const element = parts[0];
      const x = parseFloat(parts[1]);
      const y = parseFloat(parts[2]);
      const z = parseFloat(parts[3]);

      if (isNaN(x) || isNaN(y) || isNaN(z)) {
        return {
          isValid: false,
          error: `Line ${atomStartIndex + i + 1}: X, Y, Z coordinates must be numbers`,
        };
      }

      // Validate element symbol
      if (!isValidElement(element)) {
        return {
          isValid: false,
          error: `Line ${atomStartIndex + i + 1}: "${element}" is not a valid element symbol`,
        };
      }

      atoms.push({ element, x, y, z });
    }

    return {
      isValid: true,
      data: {
        numAtoms,
        comment,
        atoms,
      },
    };
  } catch (error) {
    return {
      isValid: false,
      error: `Parse error: ${error instanceof Error ? error.message : String(error)}`,
    };
  }
}

/**
 * Validate if a string is a valid chemical element symbol
 */
function isValidElement(element: string): boolean {
  const validElements = [
    'H',
    'He',
    'Li',
    'Be',
    'B',
    'C',
    'N',
    'O',
    'F',
    'Ne',
    'Na',
    'Mg',
    'Al',
    'Si',
    'P',
    'S',
    'Cl',
    'Ar',
    'K',
    'Ca',
    'Sc',
    'Ti',
    'V',
    'Cr',
    'Mn',
    'Fe',
    'Co',
    'Ni',
    'Cu',
    'Zn',
    'Ga',
    'Ge',
    'As',
    'Se',
    'Br',
    'Kr',
    'Rb',
    'Sr',
    'Y',
    'Zr',
    'Nb',
    'Mo',
    'Tc',
    'Ru',
    'Rh',
    'Pd',
    'Ag',
    'Cd',
    'In',
    'Sn',
    'Sb',
    'Te',
    'I',
    'Xe',
    'Cs',
    'Ba',
    'La',
    'Ce',
    'Pr',
    'Nd',
    'Pm',
    'Sm',
    'Eu',
    'Gd',
    'Tb',
    'Dy',
    'Ho',
    'Er',
    'Tm',
    'Yb',
    'Lu',
    'Hf',
    'Ta',
    'W',
    'Re',
    'Os',
    'Ir',
    'Pt',
    'Au',
    'Hg',
    'Tl',
    'Pb',
    'Bi',
    'Po',
    'At',
    'Rn',
    'Fr',
    'Ra',
    'Ac',
    'Th',
    'Pa',
    'U',
    'Np',
    'Pu',
    'Am',
    'Cm',
    'Bk',
    'Cf',
    'Es',
    'Fm',
    'Md',
    'No',
    'Lr',
    'Rf',
    'Db',
    'Sg',
    'Bh',
    'Hs',
    'Mt',
    'Ds',
    'Rg',
    'Cn',
    'Nh',
    'Fl',
    'Mc',
    'Lv',
    'Ts',
    'Og',
  ];

  return validElements.includes(element);
}
