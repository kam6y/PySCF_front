import {
  useRef,
  useState,
  useEffect,
  useCallback,
  useMemo,
  type ChangeEvent,
  type KeyboardEvent,
} from 'react';
import { StyleSpec, ExtendedStyleSpec } from '../../types/3dmol';
import { searchPubChem, convertSmilesToXyz } from '../api/molecule';
import { useSupportedParameters } from './useCalculationQueries';
import { useAppSettings } from './useAppSettings';
import { useGpu4PyscfStatus } from './useGpu4Pyscf';
import { useMethodDefaults } from './useMethodDefaults';
import {
  QuantumCalculationRequest,
  CalculationInstance,
  PubChemSearchResponseData,
  SMILESConvertResponseData,
} from '../types/api-types';
import { isCustomDielectricConstant } from '../constants/calculationDefaults';
import { showErrorNotification } from '../store/notificationStore';

export type DistributiveKeyOf<T> = T extends any ? keyof T : never;

interface RenameModalState {
  isOpen: boolean;
  pendingName: string;
  targetCalculationId: string;
}

export interface UseCalculationFormProps {
  activeCalculation?: CalculationInstance;
  onCalculationUpdate: (updatedCalculation: CalculationInstance) => void;
  onStartCalculation: (
    params: QuantumCalculationRequest
  ) => Promise<CalculationInstance>;
  onCalculationRename: (id: string, newName: string) => Promise<void>;
  createNewCalculationFromExisting: (
    originalCalc: CalculationInstance,
    newParams: QuantumCalculationRequest
  ) => void;
}

export const useCalculationForm = ({
  activeCalculation,
  onCalculationUpdate,
  onStartCalculation,
  onCalculationRename,
  createNewCalculationFromExisting,
}: UseCalculationFormProps) => {
  const previousCalculationIdRef = useRef<string | null>(null);
  const [currentStyle, setCurrentStyle] = useState<ExtendedStyleSpec | null>(
    null
  );
  const [calculationError, setCalculationError] = useState<string | null>(null);
  const [inputMethod, setInputMethod] = useState('pubchem');
  const [pubchemInput, setPubchemInput] = useState('');
  const [isConverting, setIsConverting] = useState(false);
  const [convertError, setConvertError] = useState<string | null>(null);
  const [localName, setLocalName] = useState('');
  const [isEditingName, setIsEditingName] = useState(false);
  const [showAxes, setShowAxes] = useState(false);
  const [showCoordinates, setShowCoordinates] = useState(false);
  const [useAtomicRadii, setUseAtomicRadii] = useState(false);
  const [isRenaming, setIsRenaming] = useState(false);
  const [renameModal, setRenameModal] = useState<RenameModalState>({
    isOpen: false,
    pendingName: '',
    targetCalculationId: '',
  });

  const {
    data: supportedParams,
    isLoading: isLoadingParams,
    error: paramsError,
  } = useSupportedParameters();
  const { applyMethodDefaults, isParameterDisabled, getParameterConstraint } =
    useMethodDefaults();
  const { settings: appSettings } = useAppSettings();
  const { data: gpuStatus } = useGpu4PyscfStatus();

  const isGpuAccelerationEnabled =
    appSettings?.gpu_acceleration_enabled ?? false;
  const gpuCapableMethods = useMemo(() => new Set(['HF', 'DFT', 'TDDFT']), []);

  const isCompleted = useMemo(
    () =>
      activeCalculation?.status === 'completed' ||
      activeCalculation?.status === 'error',
    [activeCalculation?.status]
  );

  const hasValidMolecule = useMemo(
    () =>
      !!(
        activeCalculation?.parameters?.xyz &&
        activeCalculation.parameters.xyz.trim() !== ''
      ),
    [activeCalculation?.parameters?.xyz]
  );

  const solventDisplayValue = useMemo((): string => {
    const solventValue = activeCalculation?.parameters?.solvent || '-';
    return isCustomDielectricConstant(solventValue) ? 'custom' : solventValue;
  }, [activeCalculation?.parameters?.solvent]);

  const customDielectricValue = useMemo((): string => {
    const solventValue = activeCalculation?.parameters?.solvent || '';
    return isCustomDielectricConstant(solventValue) ? solventValue : '';
  }, [activeCalculation?.parameters?.solvent]);

  const selectedMethod = activeCalculation?.parameters?.calculation_method;
  const isGpuCapableMethod = gpuCapableMethods.has(selectedMethod || 'DFT');
  const isGpuReady =
    !!gpuStatus &&
    gpuStatus.is_linux &&
    gpuStatus.cuda_detected &&
    gpuStatus.cuda_supported &&
    gpuStatus.gpu4pyscf_installed;
  const showCpuSettings =
    !isGpuAccelerationEnabled || !isGpuCapableMethod || !isGpuReady;

  useEffect(() => {
    const currentCalculationId = activeCalculation?.id || null;
    const previousCalculationId = previousCalculationIdRef.current;
    const isNewCalculation = currentCalculationId !== previousCalculationId;

    if (activeCalculation) {
      if (isNewCalculation || !isEditingName) {
        setLocalName(
          activeCalculation.name ||
            activeCalculation.parameters?.molecule_name ||
            ''
        );
      }
    } else {
      setLocalName('');
      setIsEditingName(false);
      setRenameModal({
        isOpen: false,
        pendingName: '',
        targetCalculationId: '',
      });
    }

    previousCalculationIdRef.current = currentCalculationId;
  }, [
    activeCalculation?.id,
    activeCalculation?.name,
    activeCalculation?.parameters?.molecule_name,
    isEditingName,
  ]);

  useEffect(() => {
    if (hasValidMolecule) {
      setCurrentStyle(prevStyle => {
        if (!prevStyle) return null;
        const style = { ...prevStyle };
        if (useAtomicRadii) {
          style._useAtomicRadii = true;
          style._baseAtomRadius = 0.3;
        } else {
          style._useAtomicRadii = false;
        }
        return style;
      });
    }
  }, [hasValidMolecule, useAtomicRadii]);

  const handleParamChange = useCallback(
    (
      field: DistributiveKeyOf<QuantumCalculationRequest>,
      value: string | number | boolean
    ) => {
      if (!activeCalculation) return;

      if (field === 'name') {
        if (activeCalculation.id.startsWith('new-calculation-')) {
          const stringValue = String(value);
          const updatedParams = {
            ...activeCalculation.parameters,
            [field]: stringValue,
          };
          onCalculationUpdate({
            ...activeCalculation,
            name: stringValue,
            parameters: updatedParams,
          });
        }
        return;
      }

      const currentParams = activeCalculation.parameters;
      let processedValue = value;
      if (field === 'solvent' && value === 'custom') {
        processedValue = '78.36';
      }

      let updatedParams: QuantumCalculationRequest;
      if (field === 'calculation_method') {
        updatedParams = applyMethodDefaults(
          currentParams,
          value as string
        ) as QuantumCalculationRequest;
      } else {
        updatedParams = {
          ...currentParams,
          [field]: processedValue,
        } as QuantumCalculationRequest;
      }

      if (isCompleted && field !== 'xyz') {
        createNewCalculationFromExisting(activeCalculation, updatedParams);
      } else {
        onCalculationUpdate({
          ...activeCalculation,
          parameters: updatedParams,
        });
      }
    },
    [
      activeCalculation,
      onCalculationUpdate,
      createNewCalculationFromExisting,
      applyMethodDefaults,
      isCompleted,
    ]
  );

  const handleStyleChange = useCallback((style: StyleSpec) => {
    setCurrentStyle(style);
  }, []);

  const handleXYZChange = useCallback(
    (xyzData: string, isValid: boolean) => {
      if (isValid && activeCalculation) {
        const updatedParams = {
          ...activeCalculation.parameters,
          xyz: xyzData,
        } as QuantumCalculationRequest;

        if (isCompleted) {
          createNewCalculationFromExisting(activeCalculation, updatedParams);
        } else {
          onCalculationUpdate({
            ...activeCalculation,
            parameters: updatedParams,
          });
        }
      }
    },
    [
      activeCalculation,
      onCalculationUpdate,
      createNewCalculationFromExisting,
      isCompleted,
    ]
  );

  const handleNameChange = useCallback((e: ChangeEvent<HTMLInputElement>) => {
    setLocalName(e.target.value);
    setIsEditingName(true);
  }, []);

  const handleClearName = useCallback(() => {
    setLocalName('');
  }, []);

  const handleNameBlur = useCallback(() => {
    const newName = localName.trim();

    if (!activeCalculation || !newName || activeCalculation.name === newName) {
      if (activeCalculation) setLocalName(activeCalculation.name);
      setIsEditingName(false);
      return;
    }

    if (activeCalculation.id.startsWith('new-calculation-')) {
      try {
        onCalculationUpdate({
          ...activeCalculation,
          name: newName,
          parameters: {
            ...activeCalculation.parameters,
          },
        });
      } catch (error) {
        console.error('Error updating new calculation name:', error);
        setLocalName(activeCalculation.name || '');
      }
      setIsEditingName(false);
      return;
    }

    setRenameModal({
      isOpen: true,
      pendingName: newName,
      targetCalculationId: activeCalculation.id,
    });
  }, [activeCalculation, localName, onCalculationUpdate]);

  const handleRenameConfirm = useCallback(async () => {
    if (isRenaming) return;

    if (
      !renameModal.pendingName ||
      !renameModal.targetCalculationId ||
      activeCalculation?.id !== renameModal.targetCalculationId
    ) {
      if (activeCalculation) {
        setLocalName(activeCalculation.name || '');
      }
      setRenameModal({
        isOpen: false,
        pendingName: '',
        targetCalculationId: '',
      });
      setIsEditingName(false);
      return;
    }

    setIsRenaming(true);
    try {
      await onCalculationRename(
        renameModal.targetCalculationId,
        renameModal.pendingName
      );
    } catch (error) {
      console.error('Error renaming calculation:', error);
      setLocalName(activeCalculation.name || '');
      showErrorNotification(
        'Failed to rename calculation',
        error instanceof Error ? error.message : 'Unknown error',
        renameModal.targetCalculationId
      );
    } finally {
      setRenameModal({
        isOpen: false,
        pendingName: '',
        targetCalculationId: '',
      });
      setIsEditingName(false);
      setIsRenaming(false);
    }
  }, [
    activeCalculation,
    isRenaming,
    onCalculationRename,
    renameModal.pendingName,
    renameModal.targetCalculationId,
  ]);

  const handleRenameCancel = useCallback(() => {
    if (activeCalculation) {
      setLocalName(activeCalculation.name || '');
    }
    setRenameModal({
      isOpen: false,
      pendingName: '',
      targetCalculationId: '',
    });
    setIsEditingName(false);
  }, [activeCalculation]);

  const handleNameKeyDown = useCallback(
    (e: KeyboardEvent<HTMLInputElement>) => {
      if (e.key === 'Enter') {
        e.currentTarget.blur();
      }
    },
    []
  );

  const handleStartCalculation = useCallback(async () => {
    if (
      !activeCalculation ||
      !activeCalculation.parameters?.xyz ||
      !activeCalculation.parameters.xyz.trim()
    ) {
      setCalculationError('A valid molecular structure is required.');
      return;
    }

    const moleculeName = localName.trim();
    if (!moleculeName) {
      setCalculationError('A molecule name is required.');
      return;
    }

    setCalculationError(null);

    const finalParams = {
      ...activeCalculation.parameters,
      name: moleculeName,
    } as QuantumCalculationRequest;

    try {
      const runningCalculation = await onStartCalculation(finalParams);
      onCalculationUpdate(runningCalculation);
    } catch (error) {
      setCalculationError(
        error instanceof Error ? error.message : 'An unknown error occurred.'
      );
      onCalculationUpdate({ ...activeCalculation, status: 'error' });
    }
  }, [activeCalculation, localName, onStartCalculation, onCalculationUpdate]);

  const handleXYZConvert = useCallback(async () => {
    if (!activeCalculation || !pubchemInput.trim()) return;

    setIsConverting(true);
    setConvertError(null);

    try {
      let data: PubChemSearchResponseData | SMILESConvertResponseData;
      let moleculeName = localName;

      if (inputMethod === 'smiles') {
        data = await convertSmilesToXyz(pubchemInput.trim());
        moleculeName = pubchemInput.trim();
      } else {
        const searchType = /^\d+$/.test(pubchemInput.trim()) ? 'cid' : 'name';
        data = await searchPubChem(pubchemInput.trim(), searchType);
        if ('compound_info' in data && data.compound_info?.iupac_name) {
          moleculeName = data.compound_info.iupac_name;
        }
      }

      setLocalName(moleculeName);

      const updatedParams = {
        ...activeCalculation.parameters,
        xyz: data.xyz,
        name: moleculeName,
      } as QuantumCalculationRequest;

      if (isCompleted) {
        createNewCalculationFromExisting(activeCalculation, updatedParams);
      } else {
        onCalculationUpdate({
          ...activeCalculation,
          name: moleculeName,
          parameters: updatedParams,
        });
      }
    } catch (error) {
      setConvertError(
        error instanceof Error
          ? error.message
          : 'An unknown error occurred during conversion.'
      );
    } finally {
      setIsConverting(false);
    }
  }, [
    activeCalculation,
    pubchemInput,
    inputMethod,
    localName,
    isCompleted,
    createNewCalculationFromExisting,
    onCalculationUpdate,
  ]);

  const handleInputMethodChange = useCallback((method: string) => {
    setInputMethod(method);
  }, []);

  const handlePubchemInputChange = useCallback(
    (value: string) => {
      setPubchemInput(value);
      if (convertError) {
        setConvertError(null);
      }
    },
    [convertError]
  );

  const handleConvertErrorClear = useCallback(() => {
    setConvertError(null);
  }, []);

  return {
    state: {
      currentStyle,
      calculationError,
      inputMethod,
      pubchemInput,
      isConverting,
      convertError,
      localName,
      isEditingName,
      showAxes,
      showCoordinates,
      useAtomicRadii,
      isRenaming,
      renameModal,
      params: activeCalculation?.parameters,
      calculationStatus: activeCalculation?.status,
      xyzInputValue: activeCalculation?.parameters?.xyz || '',
    },
    handlers: {
      handleParamChange,
      handleStyleChange,
      handleXYZChange,
      handleNameChange,
      handleNameBlur,
      handleNameKeyDown,
      handleStartCalculation,
      handleXYZConvert,
      handleRenameConfirm,
      handleRenameCancel,
      handleClearName,
      handleInputMethodChange,
      handlePubchemInputChange,
      handleConvertErrorClear,
      setShowAxes,
      setShowCoordinates,
      setUseAtomicRadii,
    },
    computed: {
      hasValidMolecule,
      isCompleted,
      solventDisplayValue,
      customDielectricValue,
      isGpuCapableMethod,
      showCpuSettings,
      supportedParams,
      isLoadingParams,
      paramsError,
      isParameterDisabled,
      getParameterConstraint,
    },
  };
};
