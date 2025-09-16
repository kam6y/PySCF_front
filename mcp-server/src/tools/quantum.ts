import { Tool } from '@modelcontextprotocol/sdk/types.js';
import { PySCFApiClient, QuantumCalculationRequest, PySCFApiError } from '../client.js';

export const startCalculationTool: Tool = {
  name: 'startCalculation',
  description: 'Start a quantum chemistry calculation using PySCF with various methods (DFT, HF, MP2, CCSD, TDDFT, CASCI, CASSCF)',
  inputSchema: {
    type: 'object',
    properties: {
      xyz: {
        type: 'string',
        description: 'XYZ format molecular structure data',
        minLength: 1,
      },
      name: {
        type: 'string',
        description: 'Name for this calculation',
        default: 'Unnamed Calculation',
      },
      calculation_method: {
        type: 'string',
        enum: ['DFT', 'HF', 'MP2', 'CCSD', 'CCSD_T', 'TDDFT', 'CASCI', 'CASSCF'],
        description: 'Quantum calculation method',
        default: 'DFT',
      },
      basis_function: {
        type: 'string',
        description: 'Basis set for calculation (e.g., STO-3G, 6-31G(d), cc-pVDZ)',
        default: '6-31G(d)',
      },
      exchange_correlation: {
        type: 'string',
        description: 'Exchange-correlation functional (e.g., B3LYP, PBE0, M06-2X). Required for DFT and TDDFT methods, automatically ignored for HF, MP2, CCSD, CCSD(T), CASCI, and CASSCF methods.',
        default: 'B3LYP',
      },
      charges: {
        type: 'integer',
        minimum: -10,
        maximum: 10,
        description: 'Molecular charge',
        default: 0,
      },
      spin: {
        type: 'integer',
        minimum: 0,
        maximum: 10,
        description: 'Spin (2S), number of unpaired electrons',
        default: 0,
      },
      solvent_method: {
        type: 'string',
        enum: ['none', 'ief-pcm', 'c-pcm', 'cosmo', 'ssvpe', 'ddcosmo'],
        description: 'Solvent calculation method',
        default: 'none',
      },
      solvent: {
        type: 'string',
        description: 'Solvent type or custom dielectric constant',
        default: '-',
      },
      cpu_cores: {
        type: 'integer',
        minimum: 1,
        maximum: 32,
        description: 'Number of CPU cores to use',
      },
      memory_mb: {
        type: 'integer',
        minimum: 512,
        maximum: 32768,
        description: 'Memory in MB',
      },
      tddft_nstates: {
        type: 'integer',
        minimum: 1,
        maximum: 50,
        description: 'Number of excited states to calculate (TDDFT only)',
        default: 10,
      },
      tddft_method: {
        type: 'string',
        enum: ['TDDFT', 'TDA'],
        description: 'TDDFT calculation method',
        default: 'TDDFT',
      },
      tddft_analyze_nto: {
        type: 'boolean',
        description: 'Perform Natural Transition Orbital analysis (TDDFT only)',
        default: false,
      },
      ncas: {
        type: 'integer',
        minimum: 1,
        maximum: 50,
        description: 'Number of active space orbitals (CASCI/CASSCF only)',
        default: 6,
      },
      nelecas: {
        type: 'integer',
        minimum: 1,
        maximum: 100,
        description: 'Number of active space electrons (CASCI/CASSCF only)',
        default: 8,
      },
      max_cycle_macro: {
        type: 'integer',
        minimum: 1,
        maximum: 200,
        description: 'Maximum CASSCF macro iterations (CASSCF only)',
        default: 50,
      },
      max_cycle_micro: {
        type: 'integer',
        minimum: 1,
        maximum: 20,
        description: 'Maximum CI solver micro iterations (CASCI/CASSCF)',
        default: 3,
      },
      natorb: {
        type: 'boolean',
        description: 'Transform to natural orbitals in active space (CASCI/CASSCF only)',
        default: true,
      },
      conv_tol: {
        type: 'number',
        minimum: 1e-12,
        maximum: 1e-3,
        description: 'Energy convergence tolerance (CASSCF only)',
        default: 0.0000001,
      },
      conv_tol_grad: {
        type: 'number',
        minimum: 1e-8,
        maximum: 1e-2,
        description: 'Gradient convergence tolerance (CASSCF only)',
        default: 0.0001,
      },
      optimize_geometry: {
        type: 'boolean',
        description: 'Whether to perform geometry optimization before the main calculation',
        default: true,
      },
    },
    required: ['xyz'],
  },
};

/**
 * Filter parameters based on calculation method to ensure theoretical compatibility
 */
function filterParametersForCalculationMethod(args: QuantumCalculationRequest): QuantumCalculationRequest {
  const filteredArgs = { ...args };
  const method = args.calculation_method;

  // DFT and TDDFT methods require exchange-correlation functional
  if (method === 'DFT' || method === 'TDDFT') {
    // Keep exchange_correlation parameter (already provided)
  }

  // HF, MP2, CCSD, and CCSD_T methods do not use exchange-correlation functional
  else if (method === 'HF' || method === 'MP2' || method === 'CCSD' || method === 'CCSD_T') {
    // Remove exchange_correlation to avoid theoretical incompatibility
    if ('exchange_correlation' in filteredArgs) {
      delete (filteredArgs as any).exchange_correlation;
    }
  }

  // CASCI and CASSCF methods do not use exchange-correlation functional
  else if (method === 'CASCI' || method === 'CASSCF') {
    // Remove exchange_correlation to avoid theoretical incompatibility
    if ('exchange_correlation' in filteredArgs) {
      delete (filteredArgs as any).exchange_correlation;
    }
  }

  return filteredArgs;
}

export async function handleStartCalculation(
  args: QuantumCalculationRequest,
  client: PySCFApiClient
) {
  // Filter parameters based on calculation method to avoid theoretical incompatibility
  const filteredArgs = filterParametersForCalculationMethod(args);

  try {
    const response = await client.startCalculation(filteredArgs);

    if (!response.success) {
      throw new Error(`Failed to start calculation`);
    }

    const calc = response.data.calculation;
    const params = calc.parameters;

    return {
      content: [
        {
          type: 'text',
          text: `✅ **Quantum chemistry calculation started**

**Calculation ID:** \`${calc.id}\`
**Calculation Name:** ${calc.name}
**Status:** ${calc.status}

**Calculation Parameters:**
- **Method:** ${params.calculation_method}
- **Basis Set:** ${params.basis_function}${params.exchange_correlation ? `
- **Exchange-Correlation Functional:** ${params.exchange_correlation}` : ''}
- **Charge:** ${params.charges}
- **Spin:** ${params.spin}
- **Solvent Effect:** ${params.solvent_method}${params.solvent !== '-' ? ` (${params.solvent})` : ''}
- **Geometry Optimization:** ${args.optimize_geometry !== false ? 'Enabled' : 'Disabled'}

${params.calculation_method === 'TDDFT' ? `**TDDFT Settings:**
- **Number of Excited States:** ${params.tddft_nstates}
- **Method:** ${params.tddft_method}
- **NTO Analysis:** ${params.tddft_analyze_nto ? 'Enabled' : 'Disabled'}

` : ''}${(params.calculation_method === 'CASCI' || params.calculation_method === 'CASSCF') ? `**${params.calculation_method} Settings:**
- **Active Space Orbitals:** ${params.ncas || 6}
- **Active Space Electrons:** ${params.nelecas || 6}
- **Natural Orbital Transformation:** ${params.natorb !== false ? 'Enabled' : 'Disabled'}
${params.calculation_method === 'CASSCF' ? `- **Max Macro Iterations:** ${params.max_cycle_macro || 50}
- **Energy Convergence Tolerance:** ${params.conv_tol || '1e-6'}
- **Gradient Convergence Tolerance:** ${params.conv_tol_grad || '1e-4'}` : ''}

` : ''}**Resources:**
- **CPU Cores:** ${params.cpu_cores || 'System Default'}
- **Memory:** ${params.memory_mb ? `${params.memory_mb} MB` : 'System Default'}

**Working Directory:** \`${calc.workingDirectory || 'N/A'}\`

The calculation is running in the background. Use \`getCalculationDetails\` to check progress.`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    
    // Enhanced error details for debugging using PySCFApiError
    let debugInfo = '';
    if (error instanceof PySCFApiError) {
      const details = error.details;
      debugInfo = `
**Debug Information:**
- HTTP Status: ${details.status || 'N/A'}
- Status Text: ${details.statusText || 'N/A'}
- URL: ${details.url || 'N/A'}
- Method: ${details.method || 'N/A'}
- Timestamp: ${details.timestamp}
- Response Data: ${JSON.stringify(details.responseData, null, 2) || 'N/A'}`;
    } else if (error instanceof Error && 'response' in error) {
      const axiosError = error as any;
      debugInfo = `
**Debug Information:**
- HTTP Status: ${axiosError.response?.status || 'N/A'}
- Response Data: ${JSON.stringify(axiosError.response?.data, null, 2) || 'N/A'}`;
    }
    
    return {
      content: [
        {
          type: 'text',
          text: `❌ Error occurred while starting calculation: ${errorMessage}

**Input Parameters:**
- Calculation Name: "${args.name || 'Unnamed Calculation'}"
- Method: ${args.calculation_method || 'DFT'}
- Basis Set: ${args.basis_function || '6-31G(d)'}${filteredArgs.exchange_correlation ? `
- Exchange-Correlation Functional: ${filteredArgs.exchange_correlation}` : ' (Exchange-correlation functional automatically filtered out for this method)'}
- Charge: ${args.charges || 0}
- Spin: ${args.spin || 0}
- CPU Cores: ${args.cpu_cores || 'auto'}
- Memory: ${args.memory_mb || 'auto'} MB${debugInfo}

**Possible Causes:**
- Invalid XYZ format or molecular structure
- **Theoretical parameter incompatibilities (automatically resolved in this version):**
  - HF, MP2, CCSD, CCSD(T): Exchange-correlation functionals not applicable (now automatically excluded)
  - CASCI, CASSCF: Exchange-correlation functionals not applicable (now automatically excluded)
  - DFT, TDDFT: Exchange-correlation functional required (kept when specified)
- **Method-specific issues:**
  - TDDFT: Number of excited states not specified
  - CASCI/CASSCF: Inappropriate active space parameters (ncas/nelecas)
  - CCSD(T): Requires significant computational resources
- Insufficient computational resources (CPU/Memory)
- Server internal error or temporary unavailability
- Conflict with existing calculations

**Solutions:**
1. **Verify molecular structure and format:**
   - Use \`validateXYZ\` to check XYZ format
   - Ensure reasonable molecular geometry
2. **Check method-specific requirements:**
   - **CCSD(T)**: Use correlation-consistent basis sets (cc-pVDZ, cc-pVTZ), ensure sufficient memory (>4GB)
   - **TDDFT**: Specify appropriate number of excited states (tddft_nstates)
   - **CASCI/CASSCF**: Set reasonable active space (ncas=6-12, nelecas=6-12 typically)
3. **System resources:**
   - Check available resources: \`getResourceStatus\`
   - Try lighter settings for testing (STO-3G basis, fewer cores)
4. **Retry and debugging:**
   - Check supported parameters: \`getSupportedParameters\`
   - Retry after existing calculations complete
   - Consider splitting large calculations into smaller steps`,
        },
      ],
      isError: true,
    };
  }
}

export const listCalculationsTool: Tool = {
  name: 'listCalculations',
  description: 'List all quantum chemistry calculations with their status and basic information',
  inputSchema: {
    type: 'object',
    properties: {},
  },
};

export async function handleListCalculations(_args: object, client: PySCFApiClient) {
  try {
    const response = await client.listCalculations();

    if (!response.success) {
      throw new Error(`Failed to retrieve calculation list`);
    }

    const { calculations, count, base_directory } = response.data;

    if (count === 0) {
      return {
        content: [
          {
            type: 'text',
            text: `📝 **Calculation History**

No calculations are currently saved.

Use \`startCalculation\` to start a new calculation.`,
          },
        ],
      };
    }

    const statusEmoji = (status: string) => {
      switch (status) {
        case 'completed': return '✅';
        case 'running': return '🔄';
        case 'error': return '❌';
        case 'pending': return '⏳';
        case 'waiting': return '⏸️';
        default: return '❓';
      }
    };

    const calculationsList = calculations
      .map((calc, index) => 
        `${index + 1}. ${statusEmoji(calc.status)} **${calc.name}** (ID: \`${calc.id}\`)
   - Status: ${calc.status}
   - Date: ${calc.date}
   - Checkpoint: ${calc.has_checkpoint ? 'Available' : 'None'}
   - Path: \`${calc.path}\``
      )
      .join('\n\n');

    return {
      content: [
        {
          type: 'text',
          text: `📝 **Calculation History** (${count} items)

**Base Directory:** \`${base_directory}\`

${calculationsList}

Use \`getCalculationDetails\` to view details for a specific calculation.`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ Error occurred while retrieving calculation list: ${errorMessage}`,
        },
      ],
      isError: true,
    };
  }
}

export const getCalculationDetailsTool: Tool = {
  name: 'getCalculationDetails',
  description: 'Get detailed information about a specific quantum chemistry calculation including results',
  inputSchema: {
    type: 'object',
    properties: {
      calculationId: {
        type: 'string',
        description: 'Unique calculation ID',
      },
    },
    required: ['calculationId'],
  },
};

export async function handleGetCalculationDetails(
  args: { calculationId: string },
  client: PySCFApiClient
) {
  try {
    const response = await client.getCalculationDetails(args.calculationId);

    if (!response.success) {
      throw new Error(`Failed to retrieve calculation details`);
    }

    const { calculation, files } = response.data;
    const params = calculation.parameters;
    const results = calculation.results;

    let statusText = '';
    switch (calculation.status) {
      case 'completed':
        statusText = '✅ Completed';
        break;
      case 'running':
        statusText = '🔄 Running';
        break;
      case 'error':
        statusText = `❌ Error${calculation.errorMessage ? `: ${calculation.errorMessage}` : ''}`;
        break;
      case 'pending':
        statusText = '⏳ Pending';
        break;
      case 'waiting':
        statusText = `⏸️ Waiting${calculation.waitingReason ? `: ${calculation.waitingReason}` : ''}`;
        break;
      default:
        statusText = `❓ ${calculation.status}`;
    }

    let resultText = '';
    if (results && calculation.status === 'completed') {
      resultText = `

**Calculation Results:**
- **SCF Energy:** ${results.scf_energy?.toFixed(8)} Hartree
- **Convergence:** ${results.converged ? 'Success' : 'Failed'}
- **Number of Atoms:** ${results.atom_count}
- **Occupied Orbitals:** ${results.num_occupied_orbitals}
- **Virtual Orbitals:** ${results.num_virtual_orbitals}
- **HOMO Index:** ${results.homo_index}
- **LUMO Index:** ${results.lumo_index}

**File Status:**
- **Checkpoint:** ${files.checkpoint_exists ? '✅' : '❌'}
- **Parameters File:** ${files.parameters_file_exists ? '✅' : '❌'}
- **Results File:** ${files.results_file_exists ? '✅' : '❌'}`;

      // Add TDDFT results if available
      if (results.excitation_energies && results.excitation_energies.length > 0) {
        resultText += `

**TDDFT Excited States (Top 5):**
${results.excitation_energies.slice(0, 5)
  .map((energy, i) => {
    const wavelength = results.excitation_wavelengths?.[i];
    const oscillator = results.oscillator_strengths?.[i];
    return `${i + 1}. ${energy.toFixed(4)} eV${wavelength ? ` (${wavelength.toFixed(1)} nm)` : ''}${oscillator ? `, f=${oscillator.toFixed(4)}` : ''}`;
  })
  .join('\n')}`;
      }

      // Add frequency analysis if available
      if (results.frequency_analysis_performed && results.vibrational_frequencies) {
        const freqCount = results.vibrational_frequencies.length;
        const imagCount = results.imaginary_frequencies_count || 0;
        resultText += `

**Vibrational Analysis:**
- **Real Frequencies:** ${freqCount} modes
- **Imaginary Frequencies:** ${imagCount} modes ${imagCount === 0 ? '(Optimized structure)' : '(Requires optimization)'}
- **Zero-Point Energy:** ${results.zero_point_energy?.toFixed(6)} Hartree
- **Gibbs Free Energy (298K):** ${results.gibbs_free_energy_298K?.toFixed(6)} Hartree`;
      }

      // Add CASCI/CASSCF results if available
      if (params.calculation_method === 'CASCI' || params.calculation_method === 'CASSCF') {
        if (results.casci_energy) {
          resultText += `

**CASCI Results:**
- **CASCI Energy:** ${results.casci_energy.toFixed(8)} Hartree`;
        }
        
        if (results.casscf_energy) {
          resultText += `

**CASSCF Results:**
- **CASSCF Energy:** ${results.casscf_energy.toFixed(8)} Hartree`;
          if (results.macro_iterations) {
            resultText += `
- **Macro Iterations:** ${results.macro_iterations}`;
          }
        }
        
        if (results.correlation_energy) {
          resultText += `
- **Correlation Energy:** ${results.correlation_energy.toFixed(8)} Hartree`;
        }

        // Natural orbital analysis
        if (results.natural_orbital_analysis?.enabled) {
          const noa = results.natural_orbital_analysis;
          resultText += `

**Natural Orbital Analysis:**
- **Strongly Occupied Orbitals:** ${noa.strongly_occupied_count} orbitals
- **Weakly Occupied Orbitals:** ${noa.weakly_occupied_count} orbitals
- **Virtual Orbitals:** ${noa.virtual_count} orbitals
- **Effective Electron Pairs:** ${noa.effective_electron_pairs?.toFixed(2) || 'N/A'}
- **Effective Unpaired Electrons:** ${noa.effective_unpaired_electrons?.toFixed(2) || 'N/A'}`;
        }

        // CI coefficient analysis
        if (results.ci_coefficient_analysis?.available) {
          const cia = results.ci_coefficient_analysis;
          resultText += `

**CI Coefficient Analysis:**
- **Leading Configuration:** ${cia.leading_contribution_percent?.toFixed(1) || 'N/A'}%
- **Multiconfigurational Character:** ${cia.multiconfigurational_character?.toFixed(1) || 'N/A'}%
- **Number of Major Configurations:** ${cia.major_configurations?.length || 'N/A'}`;
        }

        // Mulliken spin analysis (for open-shell systems)
        if (results.mulliken_spin_analysis?.available) {
          const msa = results.mulliken_spin_analysis;
          resultText += `

**Mulliken Spin Analysis:**
- **Total Spin Density:** ${msa.total_spin_density?.toFixed(4) || 'N/A'}
- **Expected Spin:** ${msa.expected_spin || 'N/A'}`;
        }
      }
    }

    return {
      content: [
        {
          type: 'text',
          text: `📊 **Calculation Details: ${calculation.name}**

**Basic Information:**
- **ID:** \`${calculation.id}\`
- **Status:** ${statusText}
- **Created:** ${new Date(calculation.createdAt).toLocaleString()}
- **Updated:** ${new Date(calculation.updatedAt).toLocaleString()}

**Calculation Parameters:**
- **Method:** ${params.calculation_method}
- **Basis Set:** ${params.basis_function}${params.exchange_correlation ? `
- **Exchange-Correlation Functional:** ${params.exchange_correlation}` : ''}
- **Charge:** ${params.charges}
- **Spin:** ${params.spin}
- **Solvent Effect:** ${params.solvent_method}${params.solvent !== '-' ? ` (${params.solvent})` : ''}

**Working Directory:** \`${calculation.workingDirectory || 'N/A'}\`${resultText}

${calculation.status === 'completed' ? '**Available Analysis:**\n- Molecular orbital information: `getOrbitals`\n- Orbital visualization: `getOrbitalCube`\n- Vibrational spectrum: `getIRSpectrum`' : ''}`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ Error occurred while retrieving calculation details: ${errorMessage}

Please verify the calculation ID. Available calculations can be checked with \`listCalculations\`.`,
        },
      ],
      isError: true,
    };
  }
}

export const getOrbitalsTool: Tool = {
  name: 'getOrbitals',
  description: 'Get molecular orbital information including energies and types for a completed calculation',
  inputSchema: {
    type: 'object',
    properties: {
      calculationId: {
        type: 'string',
        description: 'Unique calculation ID',
      },
    },
    required: ['calculationId'],
  },
};

export async function handleGetOrbitals(args: { calculationId: string }, client: PySCFApiClient) {
  try {
    const response = await client.getOrbitals(args.calculationId);

    if (!response.success) {
      throw new Error(`Failed to retrieve orbital information`);
    }

    const { orbitals, homo_index, lumo_index, total_orbitals, num_occupied, num_virtual } = response.data;

    // Get key orbitals around HOMO-LUMO gap
    const keyOrbitals = orbitals.filter(orbital => {
      const relativeIndex = orbital.index - homo_index;
      return relativeIndex >= -3 && relativeIndex <= 3; // HOMO-3 to LUMO+3
    });

    const orbitalList = keyOrbitals
      .map(orbital => {
        const energy_ev = orbital.energy_ev.toFixed(4);
        const occupancy = orbital.occupancy.toFixed(1);
        const typeIcon = orbital.orbital_type === 'homo' ? '🔴' : 
                        orbital.orbital_type === 'lumo' ? '🔵' : 
                        orbital.occupancy > 0 ? '⚫' : '⚪';
        return `${typeIcon} ${orbital.label}: ${energy_ev} eV (Occupancy: ${occupancy})`;
      })
      .join('\n');

    const gapEnergy = (orbitals[lumo_index]?.energy_ev || 0) - (orbitals[homo_index]?.energy_ev || 0);

    return {
      content: [
        {
          type: 'text',
          text: `🎯 **Molecular Orbital Information**

**Calculation ID:** \`${args.calculationId}\`

**Orbital Statistics:**
- **Total Orbitals:** ${total_orbitals}
- **Occupied Orbitals:** ${num_occupied}
- **Virtual Orbitals:** ${num_virtual}
- **HOMO Index:** ${homo_index}
- **LUMO Index:** ${lumo_index}
- **HOMO-LUMO Gap:** ${gapEnergy.toFixed(4)} eV (${(gapEnergy * 27.211).toFixed(1)} kcal/mol)

**Key Orbitals (around HOMO):**
${orbitalList}

**Orbital Visualization:**
Use \`getOrbitalCube\` to generate CUBE files for specific orbitals.
Example: Use \`getOrbitalCube\` to visualize HOMO (orbital ${homo_index}) or LUMO (orbital ${lumo_index})`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ Error occurred while retrieving orbital information: ${errorMessage}

**Possible Causes:**
- Calculation not completed
- Orbital data not available
- Invalid calculation ID

Check the calculation status with \`getCalculationDetails\`.`,
        },
      ],
      isError: true,
    };
  }
}

export const getOrbitalCubeTool: Tool = {
  name: 'getOrbitalCube',
  description: 'Generate CUBE file data for molecular orbital visualization',
  inputSchema: {
    type: 'object',
    properties: {
      calculationId: {
        type: 'string',
        description: 'Unique calculation ID',
      },
      orbitalIndex: {
        type: 'integer',
        minimum: 0,
        description: 'Molecular orbital index to visualize',
      },
      gridSize: {
        type: 'integer',
        minimum: 40,
        maximum: 120,
        description: 'Grid size for CUBE file generation (default 80)',
        default: 80,
      },
      isovaluePos: {
        type: 'number',
        minimum: 0.001,
        maximum: 0.1,
        description: 'Positive isovalue for visualization (default 0.02)',
        default: 0.02,
      },
      isovalueNeg: {
        type: 'number',
        minimum: -0.1,
        maximum: -0.001,
        description: 'Negative isovalue for visualization (default -0.02)',
        default: -0.02,
      },
    },
    required: ['calculationId', 'orbitalIndex'],
  },
};

export async function handleGetOrbitalCube(
  args: {
    calculationId: string;
    orbitalIndex: number;
    gridSize?: number;
    isovaluePos?: number;
    isovalueNeg?: number;
  },
  client: PySCFApiClient
) {
  try {
    const options: any = {};
    if (args.gridSize !== undefined) options.gridSize = args.gridSize;
    if (args.isovaluePos !== undefined) options.isovaluePos = args.isovaluePos;
    if (args.isovalueNeg !== undefined) options.isovalueNeg = args.isovalueNeg;
    
    const response = await client.getOrbitalCube(args.calculationId, args.orbitalIndex, options);

    if (!response.success) {
      throw new Error(`Failed to generate CUBE file`);
    }

    const { cube_data, orbital_info, generation_params, file_path, cached } = response.data;
    const cubeLines = cube_data.split('\n').length;
    const fileSizeKB = (generation_params.file_size_kb || 0).toFixed(1);

    return {
      content: [
        {
          type: 'text',
          text: `📦 **Orbital CUBE File Generation Complete**

**Orbital Information:**
- **Orbital:** ${orbital_info.label} (Index: ${orbital_info.index})
- **Energy:** ${orbital_info.energy_ev.toFixed(4)} eV
- **Occupancy:** ${orbital_info.occupancy}
- **Orbital Type:** ${orbital_info.orbital_type}

**Generation Parameters:**
- **Grid Size:** ${generation_params.grid_size}³
- **Positive Isovalue:** ${generation_params.isovalue_positive}
- **Negative Isovalue:** ${generation_params.isovalue_negative}
- **File Size:** ${fileSizeKB} KB
- **Data Lines:** ${cubeLines} lines

**File Status:**
- **Save Path:** ${file_path || 'In Memory'}
- **Cache:** ${cached ? 'Loaded from cache' : 'Newly generated'}

**CUBE Data (First 200 characters):**
\`\`\`
${cube_data.substring(0, 200)}...
\`\`\`

This CUBE file can be opened with molecular visualization software (VMD, ChemCraft, Gaussian, etc.).`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ Error occurred while generating CUBE file: ${errorMessage}

**Possible Causes:**
- Invalid orbital index
- Calculation not completed
- Checkpoint file not found

Check available orbitals with \`getOrbitals\`.`,
        },
      ],
      isError: true,
    };
  }
}

export const getIRSpectrumTool: Tool = {
  name: 'getIRSpectrum',
  description: 'Generate theoretical IR spectrum from vibrational frequency calculation',
  inputSchema: {
    type: 'object',
    properties: {
      calculationId: {
        type: 'string',
        description: 'Unique calculation ID',
      },
      broadening_fwhm: {
        type: 'number',
        minimum: 0.1,
        maximum: 1000.0,
        description: 'Full width at half maximum for Lorentzian broadening in cm⁻¹ (default 100.0)',
        default: 100.0,
      },
      x_min: {
        type: 'number',
        minimum: 0.0,
        description: 'Minimum wavenumber for spectrum range in cm⁻¹ (default 400.0)',
        default: 400.0,
      },
      x_max: {
        type: 'number',
        maximum: 10000.0,
        description: 'Maximum wavenumber for spectrum range in cm⁻¹ (default 4000.0)',
        default: 4000.0,
      },
      show_peaks: {
        type: 'boolean',
        description: 'Whether to mark individual peaks in the plot (default true)',
        default: true,
      },
    },
    required: ['calculationId'],
  },
};

export async function handleGetIRSpectrum(
  args: {
    calculationId: string;
    broadening_fwhm?: number;
    x_min?: number;
    x_max?: number;
    show_peaks?: boolean;
  },
  client: PySCFApiClient
) {
  try {
    const options: any = {};
    if (args.broadening_fwhm !== undefined) options.broadening_fwhm = args.broadening_fwhm;
    if (args.x_min !== undefined) options.x_min = args.x_min;
    if (args.x_max !== undefined) options.x_max = args.x_max;
    if (args.show_peaks !== undefined) options.show_peaks = args.show_peaks;
    
    const response = await client.getIRSpectrum(args.calculationId, options);

    if (!response.success) {
      throw new Error(`Failed to generate IR spectrum`);
    }

    const { spectrum, generation_info } = response.data;
    const metadata = spectrum.metadata;

    // Get main peaks (top 10 by intensity)
    const mainPeaks = spectrum.peaks
      .sort((a, b) => b.intensity - a.intensity)
      .slice(0, 10)
      .map((peak, i) => 
        `${(i + 1).toString().padStart(2)}. ${peak.frequency_cm.toFixed(1)} cm⁻¹ (Intensity: ${peak.intensity.toFixed(1)})`
      )
      .join('\n');

    const spectrumRange = `${args.x_min || 400} - ${args.x_max || 4000} cm⁻¹`;

    return {
      content: [
        {
          type: 'text',
          text: `📈 **Theoretical IR Spectrum**

**Calculation ID:** \`${args.calculationId}\`

**Calculation Conditions:**
- **Method:** ${metadata.method}
- **Basis Set:** ${metadata.basis_set}
- **Scale Factor:** ${metadata.scale_factor} (${metadata.scale_message})

**Spectrum Settings:**
- **Frequency Range:** ${spectrumRange}
- **Broadening (FWHM):** ${metadata.broadening_fwhm_cm} cm⁻¹
- **Number of Data Points:** ${metadata.num_points}
- **Peak Display:** ${generation_info.peaks_marked ? 'Enabled' : 'Disabled'}

**Vibrational Analysis Results:**
- **Total Frequencies:** ${metadata.num_peaks_total} modes
- **Frequencies in Range:** ${metadata.num_peaks_in_range} modes

**Major Peaks (by Intensity):**
${mainPeaks}

**Spectrum Data:**
- **X-axis (cm⁻¹):** ${spectrum.x_axis.length} points
- **Y-axis (Intensity):** Maximum ${Math.max(...spectrum.y_axis).toFixed(2)}

This spectrum can be used for comparison with experimental IR spectra. Vibrational mode analysis is useful for peak assignment.`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    
    // Enhanced error details for debugging using PySCFApiError
    let debugInfo = '';
    let specificAdvice = 'Please run calculations that include geometry optimization and vibrational analysis.';
    
    if (error instanceof PySCFApiError) {
      const details = error.details;
      debugInfo = `
**Debug Information:**
- HTTP Status: ${details.status || 'N/A'}
- Status Text: ${details.statusText || 'N/A'}
- URL: ${details.url || 'N/A'}
- Method: ${details.method || 'N/A'}
- Timestamp: ${details.timestamp}
- Calculation ID: ${args.calculationId}
- Response Data: ${JSON.stringify(details.responseData, null, 2) || 'N/A'}`;

      // Status-specific advice
      if (details.status === 404) {
        specificAdvice = `Calculation ID "${args.calculationId}" not found or IR spectrum data does not exist.\n\n**Next Steps:**\n1. Use \`getCalculationDetails\` to check calculation existence and completion status\n2. Run new calculation including vibrational analysis\n3. Ensure calculation settings include geometry optimization and vibrational analysis`;
      } else if (details.status === 400) {
        specificAdvice = 'Invalid parameters or calculation not completed. Please wait for calculation completion and retry.';
      } else if (details.status === 500) {
        specificAdvice = 'Server internal error occurred. Possible calculation data corruption or processing error.';
      }
    }
    
    return {
      content: [
        {
          type: 'text',
          text: `❌ Error occurred while generating IR spectrum: ${errorMessage}

**Possible Causes:**
- Vibrational analysis not performed
- Calculation not completed
- Vibrational frequency data not available
- Invalid calculation ID

${specificAdvice}${debugInfo}`,
        },
      ],
      isError: true,
    };
  }
}

export const getOptimizedGeometryTool: Tool = {
  name: 'getOptimizedGeometry',
  description: 'Get optimized molecular geometry from a completed calculation for use in stepwise calculations',
  inputSchema: {
    type: 'object',
    properties: {
      calculationId: {
        type: 'string',
        description: 'Unique calculation ID',
      },
    },
    required: ['calculationId'],
  },
};

export async function handleGetOptimizedGeometry(
  args: { calculationId: string },
  client: PySCFApiClient
) {
  try {
    const response = await client.getCalculationDetails(args.calculationId);

    if (!response.success) {
      throw new Error(`Failed to retrieve calculation details`);
    }

    const { calculation } = response.data;
    const params = calculation.parameters;
    const results = calculation.results;

    // Check if calculation is completed
    if (calculation.status !== 'completed') {
      return {
        content: [
          {
            type: 'text',
            text: `❌ **Optimized Coordinates Retrieval Error**

Calculation not yet completed.

**Calculation ID:** \`${calculation.id}\`
**Current Status:** ${calculation.status}

Please retry after calculation completion. Progress can be checked with \`getCalculationDetails\`.`,
          },
        ],
        isError: true,
      };
    }

    // Check if optimization was performed
    const hasOptimizedGeometry = results?.optimized_geometry;
    if (!hasOptimizedGeometry) {
      return {
        content: [
          {
            type: 'text',
            text: `❌ **Optimized Coordinates Not Available**

**Calculation ID:** \`${calculation.id}\`
**Calculation Name:** ${calculation.name}

This calculation did not perform geometry optimization or optimized coordinates are not saved.

**Possible Causes:**
- Geometry optimization was disabled
- Calculation stopped with error before geometry optimization
- Only single-point calculation was performed

**Solution:**
Run a new calculation with geometry optimization enabled:
\`startCalculation(optimize_geometry=true, ...)\``,
          },
        ],
        isError: true,
      };
    }

    // Check if geometry optimization was actually performed
    const geometryOptimized = results?.frequency_analysis_performed ||
                             (results?.imaginary_frequencies_count !== undefined);
    const optimizationStatus = geometryOptimized
      ? (results?.imaginary_frequencies_count === 0 ? 'Fully Optimized' : 'Requires Additional Optimization')
      : 'Optimization Status Unknown';

    return {
      content: [
        {
          type: 'text',
          text: `✅ **Retrieved Optimized Coordinates**

**Original Calculation Information:**
- **Calculation ID:** \`${calculation.id}\`
- **Calculation Name:** ${calculation.name}
- **Method:** ${params.calculation_method}
- **Basis Set:** ${params.basis_function}
- **Optimization Status:** ${optimizationStatus}

**Optimized Coordinates (XYZ format):**
\`\`\`
${results.optimized_geometry}
\`\`\`

**Usage:**
You can start the next calculation using these coordinates:
\`\`\`
startCalculation({
  xyz: "${results.optimized_geometry?.replace(/\n/g, '\\n')}",
  calculation_method: "CASSCF",
  optimize_geometry: false,
  // Other parameters...
})
\`\`\`

**Note:** For stepwise calculations, single-point calculations (optimize_geometry=false) are typically performed with optimized coordinates.`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ Error occurred while retrieving optimized coordinates: ${errorMessage}

Please verify the calculation ID. Available calculations can be checked with \`listCalculations\`.`,
        },
      ],
      isError: true,
    };
  }
}

