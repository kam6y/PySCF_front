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

// Calculation Instance Management Types
export interface CalculationInstance {
  id: string;                      // Unique calculation ID (directory name)
  name: string;                    // Display name for the calculation
  status: 'pending' | 'running' | 'completed' | 'error'; // Current calculation status
  createdAt: string;               // ISO timestamp
  updatedAt: string;               // ISO timestamp
  parameters: CalculationParameters;   // Input parameters
  results?: CalculationResults;    // Results (if completed)
  workingDirectory?: string;       // Path to calculation files
  errorMessage?: string;           // Error details (if status is 'error')
}

// API Response for calculation list
export interface CalculationListResponse {
  success: boolean;
  data: {
    base_directory: string;
    calculations: Array<{
      id: string;
      name: string;
      path: string;
      date: string;
      has_checkpoint: boolean;
      status: 'pending' | 'running' | 'completed' | 'error';
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
      parameters_file_exists: boolean;
      results_file_exists: boolean;
    };
  };
  error?: string;
}

// API Response for rename action
export interface RenameResponse {
  message: string;
  old_id: string;
  new_id: string;
  new_name: string;
}

// API Response for delete action
export interface DeleteResponse {
  message: string;
  deleted_id: string;
}