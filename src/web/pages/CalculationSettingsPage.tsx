import styles from './CalculationSettingsPage.module.css';
import { MoleculeViewerSection } from '../components/MoleculeViewerSection';
import { ConfirmationModal } from '../components/ConfirmationModal';
import {
  useCalculationForm,
  type UseCalculationFormProps,
} from '../hooks/useCalculationForm';
import { CalculationInstance } from '../types/api-types';
import { CalculationHeader } from '../components/calculation-settings/CalculationHeader';
import { BasicSettingsSection } from '../components/calculation-settings/BasicSettingsSection';
import { AdvancedMethodSettings } from '../components/calculation-settings/AdvancedMethodSettings';
import { SolventSettingsSection } from '../components/calculation-settings/SolventSettingsSection';
import { MolecularInputSection } from '../components/calculation-settings/MolecularInputSection';

interface CalculationSettingsPageProps extends UseCalculationFormProps {
  onCalculationPause: (id: string) => Promise<void>;
  onCalculationResume: (id: string) => Promise<void>;
}

export const CalculationSettingsPage = ({
  activeCalculation,
  onCalculationUpdate,
  onStartCalculation,
  onCalculationRename,
  onCalculationPause,
  onCalculationResume,
  createNewCalculationFromExisting,
}: CalculationSettingsPageProps) => {
  const { state, handlers, computed } = useCalculationForm({
    activeCalculation,
    onCalculationUpdate,
    onStartCalculation,
    onCalculationRename,
    createNewCalculationFromExisting,
  });

  if (!activeCalculation || !state.params || !state.calculationStatus)
    return null;

  const params = state.params;
  const calculationStatus =
    state.calculationStatus as CalculationInstance['status'];
  const sharedSectionProps = {
    params,
    calculationStatus,
    isLoadingParams: computed.isLoadingParams,
    paramsError: computed.paramsError,
    supportedParams: computed.supportedParams,
    onParamChange: handlers.handleParamChange,
  };

  return (
    <>
      <div className={styles.calculationSettingsContainers}>
        <div className={styles.calculationSettingsContainer}>
          <div className={styles.calculationHeader}>
            <CalculationHeader
              activeCalculation={activeCalculation}
              localName={state.localName}
              isEditingName={state.isEditingName}
              calculationStatus={calculationStatus}
              showCpuSettings={computed.showCpuSettings}
              hasValidMolecule={computed.hasValidMolecule}
              calculationError={state.calculationError}
              params={params}
              onNameChange={handlers.handleNameChange}
              onNameBlur={handlers.handleNameBlur}
              onNameKeyDown={handlers.handleNameKeyDown}
              onClearName={handlers.handleClearName}
              onParamChange={handlers.handleParamChange}
              onStartCalculation={handlers.handleStartCalculation}
              onCalculationPause={onCalculationPause}
              onCalculationResume={onCalculationResume}
            />
            <div className={styles.calculationColumn}>
              <BasicSettingsSection {...sharedSectionProps} />
              <AdvancedMethodSettings
                {...sharedSectionProps}
                isParameterDisabled={computed.isParameterDisabled}
              />
              <SolventSettingsSection
                {...sharedSectionProps}
                solventDisplayValue={computed.solventDisplayValue}
                customDielectricValue={computed.customDielectricValue}
              />
            </div>
          </div>
          <MoleculeViewerSection
            hasValidMolecule={computed.hasValidMolecule}
            xyzData={params.xyz}
            currentStyle={state.currentStyle}
            onStyleChange={handlers.handleStyleChange}
            showAxes={state.showAxes}
            onShowAxesChange={handlers.setShowAxes}
            showCoordinates={state.showCoordinates}
            onShowCoordinatesChange={handlers.setShowCoordinates}
            useAtomicRadii={state.useAtomicRadii}
            onUseAtomicRadiiChange={handlers.setUseAtomicRadii}
          />
        </div>
        <MolecularInputSection
          inputMethod={state.inputMethod}
          pubchemInput={state.pubchemInput}
          isConverting={state.isConverting}
          convertError={state.convertError}
          calculationStatus={calculationStatus}
          xyzInputValue={state.xyzInputValue}
          onInputMethodChange={handlers.handleInputMethodChange}
          onPubchemInputChange={handlers.handlePubchemInputChange}
          onXYZConvert={handlers.handleXYZConvert}
          onXYZChange={handlers.handleXYZChange}
          onConvertErrorClear={handlers.handleConvertErrorClear}
        />
      </div>
      <ConfirmationModal
        isOpen={state.renameModal.isOpen}
        title="Rename Calculation"
        message={`Are you sure you want to rename this calculation to "${state.renameModal.pendingName}"?`}
        confirmButtonText="Rename"
        cancelButtonText="Cancel"
        isLoading={state.isRenaming}
        variant="default"
        onConfirm={handlers.handleRenameConfirm}
        onCancel={handlers.handleRenameCancel}
      />
    </>
  );
};
