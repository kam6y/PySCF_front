// Angles in degrees. Degenerate inputs return 0 (not NaN).
// Dihedral uses the atan2 formulation from Blondel & Karplus (1996) for
// numerical stability near planar geometries.

export type Vec3 = { x: number; y: number; z: number };

const EPSILON = 1e-12;

const sub = (a: Vec3, b: Vec3): Vec3 => ({
  x: a.x - b.x,
  y: a.y - b.y,
  z: a.z - b.z,
});

const add = (a: Vec3, b: Vec3): Vec3 => ({
  x: a.x + b.x,
  y: a.y + b.y,
  z: a.z + b.z,
});

const scale = (v: Vec3, s: number): Vec3 => ({
  x: v.x * s,
  y: v.y * s,
  z: v.z * s,
});

const dot = (a: Vec3, b: Vec3): number => a.x * b.x + a.y * b.y + a.z * b.z;

const cross = (a: Vec3, b: Vec3): Vec3 => ({
  x: a.y * b.z - a.z * b.y,
  y: a.z * b.x - a.x * b.z,
  z: a.x * b.y - a.y * b.x,
});

const norm = (v: Vec3): number => Math.sqrt(dot(v, v));

const normalize = (v: Vec3): Vec3 => {
  const n = norm(v);
  if (n < EPSILON) return { x: 0, y: 0, z: 0 };
  return { x: v.x / n, y: v.y / n, z: v.z / n };
};

export const distance = (a: Vec3, b: Vec3): number => norm(sub(a, b));

export const angle = (a: Vec3, b: Vec3, c: Vec3): number => {
  const ba = sub(a, b);
  const bc = sub(c, b);
  const nBA = norm(ba);
  const nBC = norm(bc);
  if (nBA < EPSILON || nBC < EPSILON) return 0;
  const sinLen = norm(cross(ba, bc));
  const cosLen = dot(ba, bc);
  return (Math.atan2(sinLen, cosLen) * 180) / Math.PI;
};

export const dihedral = (a: Vec3, b: Vec3, c: Vec3, d: Vec3): number => {
  const b1 = sub(b, a);
  const b2 = sub(c, b);
  const b3 = sub(d, c);
  const nB2 = norm(b2);
  if (norm(b1) < EPSILON || nB2 < EPSILON || norm(b3) < EPSILON) return 0;
  const b2n = scale(b2, 1 / nB2);
  const n1 = cross(b1, b2);
  const n2 = cross(b2, b3);
  if (norm(n1) < EPSILON || norm(n2) < EPSILON) return 0;
  const m1 = cross(n1, b2n);
  const x = dot(n1, n2);
  const y = dot(m1, n2);
  return (Math.atan2(y, x) * 180) / Math.PI;
};

export const midpoint = (a: Vec3, b: Vec3): Vec3 => scale(add(a, b), 0.5);

// Falls back to the plain midpoint for collinear geometries so the label
// stays placed instead of jumping to NaN.
export const angleLabelPoint = (
  a: Vec3,
  b: Vec3,
  c: Vec3,
  offset: number
): Vec3 => {
  const mid = midpoint(a, c);
  const ba = sub(a, b);
  const bc = sub(c, b);
  const n = cross(ba, bc);
  if (norm(n) < EPSILON) return mid;
  return add(mid, scale(normalize(n), offset));
};

// Fallbacks: second-plane normal, then plain midpoint, so collinear
// geometries don't emit NaN label positions.
export const dihedralLabelPoint = (
  a: Vec3,
  b: Vec3,
  c: Vec3,
  d: Vec3,
  offset: number
): Vec3 => {
  const mid = midpoint(a, d);
  const b1 = sub(b, a);
  const b2 = sub(c, b);
  let n = cross(b1, b2);
  if (norm(n) < EPSILON) {
    const b3 = sub(d, c);
    n = cross(b2, b3);
  }
  if (norm(n) < EPSILON) return mid;
  return add(mid, scale(normalize(n), offset));
};
