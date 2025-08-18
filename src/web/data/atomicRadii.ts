/**
 * Van der Waals atomic radii in Angstrom (Å)
 * Data sourced from authoritative chemical references
 */

export interface AtomicRadiusData {
  symbol: string;
  name: string;
  vdwRadius: number; // Van der Waals radius in Å
}

// Van der Waals radii database (in Angstrom)
export const VAN_DER_WAALS_RADII: Record<string, number> = {
  // Period 1
  H: 1.2, // Hydrogen
  He: 1.4, // Helium

  // Period 2
  Li: 1.82, // Lithium
  Be: 1.53, // Beryllium
  B: 1.92, // Boron
  C: 1.7, // Carbon
  N: 1.55, // Nitrogen
  O: 1.52, // Oxygen
  F: 1.47, // Fluorine
  Ne: 1.54, // Neon

  // Period 3
  Na: 2.27, // Sodium
  Mg: 1.73, // Magnesium
  Al: 1.84, // Aluminum
  Si: 2.1, // Silicon
  P: 1.8, // Phosphorus
  S: 1.8, // Sulfur
  Cl: 1.75, // Chlorine
  Ar: 1.88, // Argon

  // Period 4
  K: 2.75, // Potassium
  Ca: 2.31, // Calcium
  Sc: 2.11, // Scandium
  Ti: 1.87, // Titanium
  V: 1.79, // Vanadium
  Cr: 1.89, // Chromium
  Mn: 1.97, // Manganese
  Fe: 1.94, // Iron
  Co: 1.92, // Cobalt
  Ni: 1.84, // Nickel
  Cu: 1.32, // Copper
  Zn: 1.22, // Zinc
  Ga: 1.87, // Gallium
  Ge: 2.11, // Germanium
  As: 1.85, // Arsenic
  Se: 1.9, // Selenium
  Br: 1.85, // Bromine
  Kr: 2.02, // Krypton

  // Period 5
  Rb: 3.03, // Rubidium
  Sr: 2.49, // Strontium
  Y: 2.32, // Yttrium
  Zr: 2.23, // Zirconium
  Nb: 2.18, // Niobium
  Mo: 2.17, // Molybdenum
  Tc: 2.16, // Technetium
  Ru: 2.13, // Ruthenium
  Rh: 2.1, // Rhodium
  Pd: 2.1, // Palladium
  Ag: 1.72, // Silver
  Cd: 1.58, // Cadmium
  In: 1.93, // Indium
  Sn: 2.17, // Tin
  Sb: 2.06, // Antimony
  Te: 2.06, // Tellurium
  I: 1.98, // Iodine
  Xe: 2.16, // Xenon

  // Period 6
  Cs: 3.43, // Cesium
  Ba: 2.68, // Barium
  La: 2.43, // Lanthanum
  Ce: 2.42, // Cerium
  Pr: 2.4, // Praseodymium
  Nd: 2.39, // Neodymium
  Pm: 2.38, // Promethium
  Sm: 2.36, // Samarium
  Eu: 2.35, // Europium
  Gd: 2.34, // Gadolinium
  Tb: 2.33, // Terbium
  Dy: 2.31, // Dysprosium
  Ho: 2.3, // Holmium
  Er: 2.29, // Erbium
  Tm: 2.27, // Thulium
  Yb: 2.26, // Ytterbium
  Lu: 2.24, // Lutetium
  Hf: 2.23, // Hafnium
  Ta: 2.22, // Tantalum
  W: 2.18, // Tungsten
  Re: 2.16, // Rhenium
  Os: 2.16, // Osmium
  Ir: 2.13, // Iridium
  Pt: 2.13, // Platinum
  Au: 1.66, // Gold
  Hg: 1.55, // Mercury
  Tl: 1.96, // Thallium
  Pb: 2.02, // Lead
  Bi: 2.07, // Bismuth
  Po: 1.97, // Polonium
  At: 2.02, // Astatine
  Rn: 2.2, // Radon

  // Period 7
  Fr: 3.48, // Francium
  Ra: 2.83, // Radium
};

/**
 * Get Van der Waals radius for an element
 */
export function getVdwRadius(element: string): number | null {
  if (!element || typeof element !== 'string' || element.trim() === '') {
    return null;
  }

  const symbol =
    element.charAt(0).toUpperCase() + element.slice(1).toLowerCase();
  return VAN_DER_WAALS_RADII[symbol] || null;
}

/**
 * Calculate relative atomic sizes based on Van der Waals radii
 * @param baseRadius - The base radius to scale from (default: hydrogen radius)
 * @param scaleFactor - Overall scaling factor for all atoms
 */
export function calculateAtomicSizes(
  baseRadius: number = VAN_DER_WAALS_RADII['H'],
  scaleFactor: number = 0.3
): Record<string, number> {
  const sizes: Record<string, number> = {};

  for (const [symbol, radius] of Object.entries(VAN_DER_WAALS_RADII)) {
    sizes[symbol] = (radius / baseRadius) * scaleFactor;
  }

  return sizes;
}

/**
 * Get atomic radius for a specific element with scaling
 */
export function getAtomicRadius(
  element: string,
  baseRadius: number = VAN_DER_WAALS_RADII['H'],
  scaleFactor: number = 0.3
): number {
  const vdwRadius = getVdwRadius(element);
  if (!vdwRadius) {
    // Fallback to default sphere radius if element not found
    return scaleFactor;
  }
  return (vdwRadius / baseRadius) * scaleFactor;
}
