import { useMemo } from 'react';
import { CalculationInstance } from '../types/api-types';

/**
 * カスタムフック: 計算結果データの処理とメモ化
 *
 * パフォーマンス最適化のため、複雑なデータ処理と条件分岐の結果をメモ化します。
 * これにより、親コンポーネントの再レンダリング時に不要な再計算を防止します。
 */
export const useProcessedCalculationResults = (
  activeCalculation?: CalculationInstance
) => {
  return useMemo(() => {
    // 計算が完了していない場合はnullを返す
    if (
      !activeCalculation ||
      activeCalculation.status !== 'completed' ||
      !activeCalculation.results
    ) {
      return null;
    }

    const results = activeCalculation.results;
    const parameters = activeCalculation.parameters;

    // 計算手法に基づいて主要なエネルギー値を決定
    let primaryEnergyLabel = '';
    let primaryEnergyValue = '';

    if (
      parameters.calculation_method === 'CCSD_T' &&
      (results as any).ccsd_t_total_energy
    ) {
      primaryEnergyLabel = 'CCSD(T) Total Energy';
      primaryEnergyValue = `${(results as any).ccsd_t_total_energy.toFixed(8)} hartree`;
    } else if (
      parameters.calculation_method === 'CCSD' &&
      (results as any).ccsd_total_energy
    ) {
      primaryEnergyLabel = 'CCSD Total Energy';
      primaryEnergyValue = `${(results as any).ccsd_total_energy.toFixed(8)} hartree`;
    } else if (
      parameters.calculation_method === 'CASSCF' &&
      (results as any).casscf_energy
    ) {
      primaryEnergyLabel = 'CASSCF Energy';
      primaryEnergyValue = `${(results as any).casscf_energy.toFixed(8)} hartree`;
    } else if (
      parameters.calculation_method === 'CASCI' &&
      (results as any).casci_energy
    ) {
      primaryEnergyLabel = 'CASCI Energy';
      primaryEnergyValue = `${(results as any).casci_energy.toFixed(8)} hartree`;
    } else {
      primaryEnergyLabel = 'SCF Energy';
      primaryEnergyValue = `${results.scf_energy?.toFixed(8) || 'N/A'} hartree`;
    }

    // 各セクションの表示判定
    const shouldShowElectronicProperties =
      (results.mulliken_charges && results.mulliken_charges.length > 0) ||
      ((results as any).mulliken_spin_analysis &&
        (results as any).mulliken_spin_analysis.available);

    const shouldShowCASSection =
      parameters.calculation_method === 'CASCI' ||
      parameters.calculation_method === 'CASSCF';

    const shouldShowTDDFTSection =
      parameters.calculation_method === 'TDDFT' && results.excitation_energies;

    const shouldShowCCSDSection =
      parameters.calculation_method === 'CCSD' ||
      parameters.calculation_method === 'CCSD_T';

    const shouldShowVibrationalSection = results.frequency_analysis_performed;

    return {
      results,
      parameters,
      primaryEnergyLabel,
      primaryEnergyValue,
      shouldShowElectronicProperties,
      shouldShowCASSection,
      shouldShowTDDFTSection,
      shouldShowCCSDSection,
      shouldShowVibrationalSection,
    };
  }, [activeCalculation]);
};
