interface NumberInputBounds {
  min: number;
  max: number;
  fallback: number;
}

export const clampIntegerInput = (
  rawValue: string,
  bounds: NumberInputBounds
): number => {
  const trimmedValue = rawValue.trim();
  if (trimmedValue === '') {
    return bounds.fallback;
  }

  const numericValue = Number(trimmedValue);
  if (!Number.isFinite(numericValue)) {
    return bounds.fallback;
  }

  const integerValue = Math.trunc(numericValue);
  return Math.min(bounds.max, Math.max(bounds.min, integerValue));
};
