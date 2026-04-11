import * as $3Dmol from '3dmol';

const SELECTION_TINT_FALLBACK = '#ffff00';
const SELECTION_TINT_DELTA = 0.4;

// delta in [-1,1]; positive lightens toward white, negative darkens toward black.
export const shadeHex = (hex: string, delta: number): string => {
  const m = hex.replace('#', '');
  if (!/^[0-9a-fA-F]{6}$/.test(m)) return SELECTION_TINT_FALLBACK;
  const r = parseInt(m.substring(0, 2), 16);
  const g = parseInt(m.substring(2, 4), 16);
  const b = parseInt(m.substring(4, 6), 16);
  const adjust = (c: number): number => {
    const target = delta >= 0 ? 255 : 0;
    return Math.round(c + (target - c) * Math.abs(delta));
  };
  const toHex = (n: number): string => n.toString(16).padStart(2, '0');
  return `#${toHex(adjust(r))}${toHex(adjust(g))}${toHex(adjust(b))}`;
};

export const elementTint = (elem: string | undefined): string => {
  if (!elem) {
    return SELECTION_TINT_FALLBACK;
  }
  const palette = $3Dmol.elementColors?.defaultColors;
  if (!palette) {
    return SELECTION_TINT_FALLBACK;
  }
  const key = elem.charAt(0).toUpperCase() + elem.slice(1).toLowerCase();
  const numColor = palette[key];
  if (numColor === undefined) {
    return SELECTION_TINT_FALLBACK;
  }
  const hex = '#' + numColor.toString(16).padStart(6, '0');
  return shadeHex(hex, SELECTION_TINT_DELTA);
};
