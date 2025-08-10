/**
 * Type definitions for quantum chemistry calculations
 */

// Calculation status types
export type CalculationStatus = 'idle' | 'running' | 'completed' | 'error';

// Calculation method types
export type CalculationMethod = 'DFT' | 'HF' | 'MP2';
export type BasisFunction = 'STO-3G' | '3-21G' | '6-31G' | '6-31G(d)' | 'cc-pVDZ';
export type ExchangeCorrelation = 'B3LYP' | 'PBE' | 'BP86' | 'wB97XD';
export type SolventMethod = 'none' | 'PCM' | 'SMD';

// Calculation parameters interface
export interface CalculationParameters {
  calculation_method: CalculationMethod;
  basis_function: BasisFunction;
  exchange_correlation: ExchangeCorrelation;
  charges: number;
  spin_multiplicity: number;
  solvent_method: SolventMethod;
  solvent: string;
  xyz: string;
  molecule_name?: string;
  cpu_cores?: number;
  memory_mb?: number;
}

// Calculation results interface
export interface CalculationResults {
  scf_energy: number;
  converged: boolean;
  homo_index: number;
  lumo_index: number;
  num_occupied_orbitals: number;
  num_virtual_orbitals: number;
  checkpoint_file: string;
  checkpoint_exists: boolean;
  working_directory: string;
  optimized_geometry: string;
  basis: string;
  xc_functional: string;
  charge: number;
  spin_multiplicity: number;
  max_cycle: number;
  atom_count: number;
}

// API response interfaces
export interface QuantumCalculationResponse {
  success: boolean;
  data: {
    calculation_results: CalculationResults;
    calculation_parameters: {
      method: CalculationMethod;
      basis: BasisFunction;
      xc_functional: ExchangeCorrelation;
      charge: number;
      spin_multiplicity: number;
    };
  };
  error?: string;
}

// Calculation state interface for component state management
export interface CalculationState {
  status: CalculationStatus;
  results: CalculationResults | null;
  parameters: CalculationParameters | null;
  error: string | null;
  startTime: Date | null;
  endTime: Date | null;
}

// Props for calculation-related components
export interface CalculationStatusProps {
  status: CalculationStatus;
  error?: string | null;
  startTime?: Date | null;
}

export interface ResultsViewerProps {
  results: CalculationResults;
  parameters: CalculationParameters;
}

export interface OrbitalDisplayProps {
  homoIndex: number;
  lumoIndex: number;
  numOccupied: number;
  numVirtual: number;
}

// Calculation Instance Management Types
export interface CalculationInstance {
  id: string;                                    // Unique calculation ID
  name: string;                                  // Display name for the calculation
  status: 'pending' | 'running' | 'completed' | 'error';  // Current calculation status
  createdAt: string;                            // ISO timestamp
  updatedAt: string;                            // ISO timestamp
  parameters: CalculationParameters;            // Input parameters
  results?: CalculationResults;                 // Results (if completed)
  workingDirectory?: string;                    // Path to calculation files
  errorMessage?: string;                        // Error details (if status is 'error')
}

// API Response for calculation list
export interface CalculationListResponse {
  success: boolean;
  data: {
    base_directory: string;
    calculations: Array<{
      name: string;
      path: string;
      date: string;
      has_checkpoint: boolean;
    }>;
    count: number;
  };
  error?: string;
}

// API Response for calculation details
export interface CalculationDetailsResponse {
  success: boolean;
  data: {
    calculation: CalculationInstance;
    files: {
      checkpoint_exists: boolean;
      info_file_exists: boolean;
      geometry_file_exists: boolean;
    };
  };
  error?: string;
}

// Calculation update request
export interface CalculationUpdateRequest {
  name?: string;
  parameters?: Partial<CalculationParameters>;
}

// Application state for calculation management
export interface CalculationManagementState {
  calculations: CalculationInstance[];
  activeCalculationId: string | null;
  isLoading: boolean;
  error: string | null;
}

// Sidebar calculation item props
export interface SidebarCalculationItemProps {
  calculation: CalculationInstance;
  isActive: boolean;
  onSelect: (id: string) => void;
  onRename: (id: string, newName: string) => void;
  onDelete: (id: string) => void;
}