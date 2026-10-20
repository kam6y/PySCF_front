import { useCallback, useEffect, useState } from 'react';

export interface UseAtomMeasurementResult {
  selectedAtomIndices: number[];
  handleAtomClick: (atomIndex: number) => void;
}

export function useAtomMeasurement(
  resetKey: string | null | undefined
): UseAtomMeasurementResult {
  const [selectedAtomIndices, setSelectedAtomIndices] = useState<number[]>([]);

  // Re-clicking a selected atom deselects it; a 5th distinct click restarts the measurement.
  const handleAtomClick = useCallback((atomIndex: number) => {
    setSelectedAtomIndices(previousIndices => {
      if (previousIndices.includes(atomIndex)) {
        return previousIndices.filter(i => i !== atomIndex);
      }
      if (previousIndices.length === 4) {
        return [atomIndex];
      }
      return [...previousIndices, atomIndex];
    });
  }, []);

  useEffect(() => {
    setSelectedAtomIndices(prev => (prev.length === 0 ? prev : []));
  }, [resetKey]);

  return {
    selectedAtomIndices,
    handleAtomClick,
  };
}
