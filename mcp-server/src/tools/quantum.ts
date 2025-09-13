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
        description: 'Exchange-correlation functional (e.g., B3LYP, PBE0, M06-2X)',
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
    },
    required: ['xyz'],
  },
};

export async function handleStartCalculation(
  args: QuantumCalculationRequest,
  client: PySCFApiClient
) {
  try {
    const response = await client.startCalculation(args);

    if (!response.success) {
      throw new Error(`計算開始に失敗しました`);
    }

    const calc = response.data.calculation;
    const params = calc.parameters;

    return {
      content: [
        {
          type: 'text',
          text: `✅ **量子化学計算を開始しました**

**計算ID:** \`${calc.id}\`
**計算名:** ${calc.name}
**ステータス:** ${calc.status}

**計算パラメータ:**
- **手法:** ${params.calculation_method}
- **基底関数:** ${params.basis_function}
- **交換相関汎関数:** ${params.exchange_correlation}
- **電荷:** ${params.charges}
- **スピン:** ${params.spin}
- **溶媒効果:** ${params.solvent_method}${params.solvent !== '-' ? ` (${params.solvent})` : ''}

${params.calculation_method === 'TDDFT' ? `**TDDFT設定:**
- **励起状態数:** ${params.tddft_nstates}
- **手法:** ${params.tddft_method}
- **NTO解析:** ${params.tddft_analyze_nto ? '有効' : '無効'}

` : ''}${(params.calculation_method === 'CASCI' || params.calculation_method === 'CASSCF') ? `**${params.calculation_method}設定:**
- **アクティブ空間軌道数:** ${params.ncas || 6}
- **アクティブ空間電子数:** ${params.nelecas || 6}
- **自然軌道変換:** ${params.natorb !== false ? '有効' : '無効'}
${params.calculation_method === 'CASSCF' ? `- **最大マクロイテレーション:** ${params.max_cycle_macro || 50}
- **エネルギー収束許容値:** ${params.conv_tol || '1e-6'}
- **勾配収束許容値:** ${params.conv_tol_grad || '1e-4'}` : ''}

` : ''}**リソース:**
- **CPUコア:** ${params.cpu_cores || 'システム設定'}
- **メモリ:** ${params.memory_mb ? `${params.memory_mb} MB` : 'システム設定'}

**作業ディレクトリ:** \`${calc.workingDirectory || 'N/A'}\`

計算は背景で実行されています。進行状況を確認するには \`getCalculationDetails\` を使用してください。`,
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
**デバッグ情報:**
- HTTPステータス: ${details.status || 'N/A'}
- ステータステキスト: ${details.statusText || 'N/A'}
- URL: ${details.url || 'N/A'}
- メソッド: ${details.method || 'N/A'}
- タイムスタンプ: ${details.timestamp}
- レスポンスデータ: ${JSON.stringify(details.responseData, null, 2) || 'N/A'}`;
    } else if (error instanceof Error && 'response' in error) {
      const axiosError = error as any;
      debugInfo = `
**デバッグ情報:**
- HTTPステータス: ${axiosError.response?.status || 'N/A'}
- レスポンスデータ: ${JSON.stringify(axiosError.response?.data, null, 2) || 'N/A'}`;
    }
    
    return {
      content: [
        {
          type: 'text',
          text: `❌ 計算開始でエラーが発生しました: ${errorMessage}

**入力パラメータ:**
- 計算名: "${args.name || 'Unnamed Calculation'}"
- 手法: ${args.calculation_method || 'DFT'}
- 基底関数: ${args.basis_function || '6-31G(d)'}
- 交換相関汎関数: ${args.exchange_correlation || 'B3LYP'}
- 電荷: ${args.charges || 0}
- スピン: ${args.spin || 0}
- CPUコア: ${args.cpu_cores || 'auto'}
- メモリ: ${args.memory_mb || 'auto'} MB${debugInfo}

**可能な原因:**
- 無効なXYZ形式
- **理論的に不適切なパラメータの組み合わせ**
  - HF法で交換相関汎関数が指定されている（HF法には交換相関汎関数は不要です）
  - TDDFT法で励起状態数が指定されていない
  - CASCI/CASSCF法でアクティブ空間パラメータが不適切
- リソース不足（CPU/メモリ）
- サーバー内部エラー
- 既存の計算との競合

**解決方法:**
1. **パラメータの理論的適合性を確認:**
   - HF法: 交換相関汎関数は不要（自動的に無視されます）
   - DFT/TDDFT法: 適切な交換相関汎関数（B3LYP, PBE0等）を指定
   - CASCI/CASSCF法: ncas（軌道数）とnelecas（電子数）を適切に設定
2. XYZ形式を\`validateXYZ\`で確認してください
3. \`getSupportedParameters\`で利用可能なパラメータを確認してください
4. \`getResourceStatus\`でシステムリソースを確認してください
5. より軽量な設定（STO-3G基底関数など）で試してください
6. 既存の計算が完了してから再試行してください`,
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
      throw new Error(`計算一覧の取得に失敗しました`);
    }

    const { calculations, count, base_directory } = response.data;

    if (count === 0) {
      return {
        content: [
          {
            type: 'text',
            text: `📝 **計算履歴**

現在保存されている計算はありません。

新しい計算を開始するには \`startCalculation\` を使用してください。`,
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
   - ステータス: ${calc.status}
   - 日付: ${calc.date}
   - チェックポイント: ${calc.has_checkpoint ? '有' : '無'}
   - パス: \`${calc.path}\``
      )
      .join('\n\n');

    return {
      content: [
        {
          type: 'text',
          text: `📝 **計算履歴** (${count}件)

**ベースディレクトリ:** \`${base_directory}\`

${calculationsList}

特定の計算の詳細を確認するには \`getCalculationDetails\` を使用してください。`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ 計算一覧の取得でエラーが発生しました: ${errorMessage}`,
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
      throw new Error(`計算詳細の取得に失敗しました`);
    }

    const { calculation, files } = response.data;
    const params = calculation.parameters;
    const results = calculation.results;

    let statusText = '';
    switch (calculation.status) {
      case 'completed':
        statusText = '✅ 完了';
        break;
      case 'running':
        statusText = '🔄 実行中';
        break;
      case 'error':
        statusText = `❌ エラー${calculation.errorMessage ? `: ${calculation.errorMessage}` : ''}`;
        break;
      case 'pending':
        statusText = '⏳ 待機中';
        break;
      case 'waiting':
        statusText = `⏸️ 待機中${calculation.waitingReason ? `: ${calculation.waitingReason}` : ''}`;
        break;
      default:
        statusText = `❓ ${calculation.status}`;
    }

    let resultText = '';
    if (results && calculation.status === 'completed') {
      resultText = `

**計算結果:**
- **SCFエネルギー:** ${results.scf_energy?.toFixed(8)} Hartree
- **収束:** ${results.converged ? '成功' : '失敗'}
- **原子数:** ${results.atom_count}
- **占有軌道数:** ${results.num_occupied_orbitals}
- **仮想軌道数:** ${results.num_virtual_orbitals}
- **HOMO指数:** ${results.homo_index}
- **LUMO指数:** ${results.lumo_index}

**ファイル状況:**
- **チェックポイント:** ${files.checkpoint_exists ? '✅' : '❌'}
- **パラメータファイル:** ${files.parameters_file_exists ? '✅' : '❌'}
- **結果ファイル:** ${files.results_file_exists ? '✅' : '❌'}`;

      // Add TDDFT results if available
      if (results.excitation_energies && results.excitation_energies.length > 0) {
        resultText += `

**TDDFT励起状態 (上位5状態):**
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

**振動解析:**
- **実振動数:** ${freqCount}個
- **虚振動数:** ${imagCount}個 ${imagCount === 0 ? '(最適化済み構造)' : '(要最適化)'}
- **零点エネルギー:** ${results.zero_point_energy?.toFixed(6)} Hartree
- **自由エネルギー (298K):** ${results.gibbs_free_energy_298K?.toFixed(6)} Hartree`;
      }

      // Add CASCI/CASSCF results if available
      if (params.calculation_method === 'CASCI' || params.calculation_method === 'CASSCF') {
        if (results.casci_energy) {
          resultText += `

**CASCI結果:**
- **CASCIエネルギー:** ${results.casci_energy.toFixed(8)} Hartree`;
        }
        
        if (results.casscf_energy) {
          resultText += `

**CASSCF結果:**
- **CASSCFエネルギー:** ${results.casscf_energy.toFixed(8)} Hartree`;
          if (results.macro_iterations) {
            resultText += `
- **マクロイテレーション数:** ${results.macro_iterations}`;
          }
        }
        
        if (results.correlation_energy) {
          resultText += `
- **相関エネルギー:** ${results.correlation_energy.toFixed(8)} Hartree`;
        }

        // Natural orbital analysis
        if (results.natural_orbital_analysis?.enabled) {
          const noa = results.natural_orbital_analysis;
          resultText += `

**自然軌道解析:**
- **強占有軌道:** ${noa.strongly_occupied_count}個
- **弱占有軌道:** ${noa.weakly_occupied_count}個  
- **仮想軌道:** ${noa.virtual_count}個
- **有効電子対数:** ${noa.effective_electron_pairs?.toFixed(2) || 'N/A'}
- **有効不対電子数:** ${noa.effective_unpaired_electrons?.toFixed(2) || 'N/A'}`;
        }

        // CI coefficient analysis
        if (results.ci_coefficient_analysis?.available) {
          const cia = results.ci_coefficient_analysis;
          resultText += `

**CI係数解析:**
- **主要配置:** ${cia.leading_contribution_percent?.toFixed(1) || 'N/A'}%
- **多配置性:** ${cia.multiconfigurational_character?.toFixed(1) || 'N/A'}%
- **有効配置数:** ${cia.major_configurations?.length || 'N/A'}`;
        }

        // Mulliken spin analysis (for open-shell systems)
        if (results.mulliken_spin_analysis?.available) {
          const msa = results.mulliken_spin_analysis;
          resultText += `

**Mulliken スピン解析:**
- **総スピン密度:** ${msa.total_spin_density?.toFixed(4) || 'N/A'}
- **期待スピン:** ${msa.expected_spin || 'N/A'}`;
        }
      }
    }

    return {
      content: [
        {
          type: 'text',
          text: `📊 **計算詳細: ${calculation.name}**

**基本情報:**
- **ID:** \`${calculation.id}\`
- **ステータス:** ${statusText}
- **作成日時:** ${new Date(calculation.createdAt).toLocaleString()}
- **更新日時:** ${new Date(calculation.updatedAt).toLocaleString()}

**計算パラメータ:**
- **手法:** ${params.calculation_method}
- **基底関数:** ${params.basis_function}
- **交換相関汎関数:** ${params.exchange_correlation}
- **電荷:** ${params.charges}
- **スピン:** ${params.spin}
- **溶媒効果:** ${params.solvent_method}${params.solvent !== '-' ? ` (${params.solvent})` : ''}

**作業ディレクトリ:** \`${calculation.workingDirectory || 'N/A'}\`${resultText}

${calculation.status === 'completed' ? '**利用可能な解析:**\n- 分子軌道情報: `getOrbitals`\n- 軌道可視化: `getOrbitalCube`\n- 振動スペクトラム: `getIRSpectrum`' : ''}`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ 計算詳細の取得でエラーが発生しました: ${errorMessage}

計算IDを確認してください。利用可能な計算は \`listCalculations\` で確認できます。`,
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
      throw new Error(`軌道情報の取得に失敗しました`);
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
        return `${typeIcon} ${orbital.label}: ${energy_ev} eV (占有: ${occupancy})`;
      })
      .join('\n');

    const gapEnergy = (orbitals[lumo_index]?.energy_ev || 0) - (orbitals[homo_index]?.energy_ev || 0);

    return {
      content: [
        {
          type: 'text',
          text: `🎯 **分子軌道情報**

**計算ID:** \`${args.calculationId}\`

**軌道統計:**
- **総軌道数:** ${total_orbitals}
- **占有軌道数:** ${num_occupied}
- **仮想軌道数:** ${num_virtual}
- **HOMO指数:** ${homo_index}
- **LUMO指数:** ${lumo_index}
- **HOMO-LUMOギャップ:** ${gapEnergy.toFixed(4)} eV (${(gapEnergy * 27.211).toFixed(1)} kcal/mol)

**主要軌道 (HOMO周辺):**
${orbitalList}

**軌道可視化:**
特定の軌道のCUBEファイルを生成するには \`getOrbitalCube\` を使用してください。
例: \`getOrbitalCube\` で HOMO (軌道${homo_index}) や LUMO (軌道${lumo_index}) を可視化`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ 軌道情報の取得でエラーが発生しました: ${errorMessage}

**可能な原因:**
- 計算が完了していない
- 軌道データが利用できない
- 計算IDが無効

計算のステータスを \`getCalculationDetails\` で確認してください。`,
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
      throw new Error(`CUBEファイル生成に失敗しました`);
    }

    const { cube_data, orbital_info, generation_params, file_path, cached } = response.data;
    const cubeLines = cube_data.split('\n').length;
    const fileSizeKB = (generation_params.file_size_kb || 0).toFixed(1);

    return {
      content: [
        {
          type: 'text',
          text: `📦 **軌道CUBEファイル生成完了**

**軌道情報:**
- **軌道:** ${orbital_info.label} (指数: ${orbital_info.index})
- **エネルギー:** ${orbital_info.energy_ev.toFixed(4)} eV
- **占有度:** ${orbital_info.occupancy}
- **軌道種別:** ${orbital_info.orbital_type}

**生成パラメータ:**
- **グリッドサイズ:** ${generation_params.grid_size}³
- **正のアイソ値:** ${generation_params.isovalue_positive}
- **負のアイソ値:** ${generation_params.isovalue_negative}
- **ファイルサイズ:** ${fileSizeKB} KB
- **データ行数:** ${cubeLines}行

**ファイル状況:**
- **保存パス:** ${file_path || 'メモリ内'}
- **キャッシュ:** ${cached ? 'キャッシュから読み込み' : '新規生成'}

**CUBEデータ (先頭200文字):**
\`\`\`
${cube_data.substring(0, 200)}...
\`\`\`

このCUBEファイルは分子可視化ソフトウェア（VMD、ChemCraft、Gaussianなど）で開くことができます。`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ CUBEファイル生成でエラーが発生しました: ${errorMessage}

**可能な原因:**
- 無効な軌道指数
- 計算が完了していない
- チェックポイントファイルが見つからない

利用可能な軌道は \`getOrbitals\` で確認してください。`,
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
      throw new Error(`IRスペクトラム生成に失敗しました`);
    }

    const { spectrum, generation_info } = response.data;
    const metadata = spectrum.metadata;

    // Get main peaks (top 10 by intensity)
    const mainPeaks = spectrum.peaks
      .sort((a, b) => b.intensity - a.intensity)
      .slice(0, 10)
      .map((peak, i) => 
        `${(i + 1).toString().padStart(2)}. ${peak.frequency_cm.toFixed(1)} cm⁻¹ (強度: ${peak.intensity.toFixed(1)})`
      )
      .join('\n');

    const spectrumRange = `${args.x_min || 400} - ${args.x_max || 4000} cm⁻¹`;

    return {
      content: [
        {
          type: 'text',
          text: `📈 **理論IRスペクトラム**

**計算ID:** \`${args.calculationId}\`

**計算条件:**
- **手法:** ${metadata.method}
- **基底関数:** ${metadata.basis_set}
- **スケール因子:** ${metadata.scale_factor} (${metadata.scale_message})

**スペクトラム設定:**
- **周波数範囲:** ${spectrumRange}
- **ブロードニング (FWHM):** ${metadata.broadening_fwhm_cm} cm⁻¹
- **データ点数:** ${metadata.num_points}
- **ピーク表示:** ${generation_info.peaks_marked ? '有効' : '無効'}

**振動解析結果:**
- **全振動数:** ${metadata.num_peaks_total}個
- **範囲内振動数:** ${metadata.num_peaks_in_range}個

**主要ピーク (強度順):**
${mainPeaks}

**スペクトラムデータ:**
- **X軸 (cm⁻¹):** ${spectrum.x_axis.length}点
- **Y軸 (強度):** 最大値 ${Math.max(...spectrum.y_axis).toFixed(2)}

このスペクトラムは実験IRスペクトラムとの比較に使用できます。ピークの帰属には振動モード解析が有効です。`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    
    // Enhanced error details for debugging using PySCFApiError
    let debugInfo = '';
    let specificAdvice = '構造最適化と振動解析を含む計算を実行してください。';
    
    if (error instanceof PySCFApiError) {
      const details = error.details;
      debugInfo = `
**デバッグ情報:**
- HTTPステータス: ${details.status || 'N/A'}
- ステータステキスト: ${details.statusText || 'N/A'}
- URL: ${details.url || 'N/A'}
- メソッド: ${details.method || 'N/A'}
- タイムスタンプ: ${details.timestamp}
- 計算ID: ${args.calculationId}
- レスポンスデータ: ${JSON.stringify(details.responseData, null, 2) || 'N/A'}`;

      // Status-specific advice
      if (details.status === 404) {
        specificAdvice = `計算ID "${args.calculationId}" が見つからないか、IRスペクトラムデータが存在しません。\n\n**次のステップ:**\n1. \`getCalculationDetails\` で計算の存在と完了状況を確認\n2. 振動解析を含む新しい計算を実行\n3. 計算設定に構造最適化と振動解析が含まれていることを確認`;
      } else if (details.status === 400) {
        specificAdvice = '無効なパラメータまたは計算が未完了です。計算の完了を待ってから再試行してください。';
      } else if (details.status === 500) {
        specificAdvice = 'サーバー内部エラーが発生しました。計算データの破損または処理エラーの可能性があります。';
      }
    }
    
    return {
      content: [
        {
          type: 'text',
          text: `❌ IRスペクトラム生成でエラーが発生しました: ${errorMessage}

**可能な原因:**
- 振動解析が実行されていない
- 計算が完了していない
- 振動数データが利用できない
- 無効な計算ID

${specificAdvice}${debugInfo}`,
        },
      ],
      isError: true,
    };
  }
}