import { Tool } from '@modelcontextprotocol/sdk/types.js';
import { PySCFApiClient, SettingsUpdateRequest, PySCFApiError } from '../client.js';

export const getSupportedParametersTool: Tool = {
  name: 'getSupportedParameters',
  description: 'Get lists of supported quantum chemistry parameters including methods, basis sets, functionals, and solvents',
  inputSchema: {
    type: 'object',
    properties: {},
  },
};

export async function handleGetSupportedParameters(_args: object, client: PySCFApiClient) {
  try {
    const response = await client.getSupportedParameters();

    if (!response.success) {
      throw new Error(`サポートパラメータの取得に失敗しました`);
    }

    const data = response.data;

    // Format basis functions by category
    const basisText = Object.entries(data.basis_functions)
      .map(([category, functions]) => `**${category}:**\n${functions.map(f => `  - ${f}`).join('\n')}`)
      .join('\n\n');

    // Format exchange-correlation functionals by category
    const xcText = Object.entries(data.exchange_correlation)
      .map(([category, functionals]) => `**${category}:**\n${functionals.map(f => `  - ${f}`).join('\n')}`)
      .join('\n\n');

    // Format solvents by category
    const solventText = Object.entries(data.solvents)
      .map(([category, solvents]) => 
        `**${category}:**\n${solvents.map(s => `  - ${s.display} (ε=${s.dielectric_constant})`).join('\n')}`
      )
      .join('\n\n');

    return {
      content: [
        {
          type: 'text',
          text: `⚙️ **サポートされる計算パラメータ**

**計算手法:**
${data.calculation_methods.map(method => `- ${method}`).join('\n')}

**基底関数:**
${basisText}

**交換相関汎関数:**
${xcText}

**溶媒効果手法:**
${data.solvent_methods.map(method => `- ${method}`).join('\n')}

**溶媒:**
${solventText}

**TDDFT手法:**
${data.tddft_methods.map(method => `- ${method}`).join('\n')}

これらのパラメータは \`startCalculation\` で使用できます。`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ サポートパラメータの取得でエラーが発生しました: ${errorMessage}

サーバーが起動していることを確認してください。`,
        },
      ],
      isError: true,
    };
  }
}

export const getSettingsTool: Tool = {
  name: 'getSettings',
  description: 'Get current application settings including parallel processing limits and resource constraints',
  inputSchema: {
    type: 'object',
    properties: {},
  },
};

export async function handleGetSettings(_args: object, client: PySCFApiClient) {
  try {
    const response = await client.getSettings();

    if (!response.success) {
      throw new Error(`設定の取得に失敗しました`);
    }

    const settings = response.data.settings;

    // Safely handle undefined settings with defaults
    const safeSettings = {
      max_parallel_instances: settings.max_parallel_instances ?? 1,
      max_cpu_utilization_percent: settings.max_cpu_utilization_percent ?? 95,
      max_memory_utilization_percent: settings.max_memory_utilization_percent ?? 95,
      system_total_cores: settings.system_total_cores ?? 1,
      system_total_memory_mb: settings.system_total_memory_mb ?? 1024
    };

    return {
      content: [
        {
          type: 'text',
          text: `⚙️ **アプリケーション設定**

**並列処理設定:**
- **最大並列インスタンス数:** ${safeSettings.max_parallel_instances}
- **最大CPU使用率:** ${safeSettings.max_cpu_utilization_percent}%
- **最大メモリ使用率:** ${safeSettings.max_memory_utilization_percent}%

**システム情報:**
- **総CPUコア数:** ${safeSettings.system_total_cores}
- **総メモリ:** ${(safeSettings.system_total_memory_mb / 1024).toFixed(1)} GB (${safeSettings.system_total_memory_mb} MB)

**実効制限:**
- **利用可能CPUコア:** ${Math.floor(safeSettings.system_total_cores * safeSettings.max_cpu_utilization_percent / 100)}コア
- **利用可能メモリ:** ${(safeSettings.system_total_memory_mb * safeSettings.max_memory_utilization_percent / 100 / 1024).toFixed(1)} GB

設定を変更するには \`updateSettings\` を使用してください。`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ 設定の取得でエラーが発生しました: ${errorMessage}`,
        },
      ],
      isError: true,
    };
  }
}

export const updateSettingsTool: Tool = {
  name: 'updateSettings',
  description: 'Update application settings for parallel processing and resource usage limits',
  inputSchema: {
    type: 'object',
    properties: {
      max_parallel_instances: {
        type: 'integer',
        minimum: 1,
        maximum: 32,
        description: 'Maximum number of parallel calculation instances',
      },
      max_cpu_utilization_percent: {
        type: 'number',
        minimum: 10.0,
        maximum: 100.0,
        description: 'Maximum CPU utilization percentage for the system',
      },
      max_memory_utilization_percent: {
        type: 'number',
        minimum: 10.0,
        maximum: 100.0,
        description: 'Maximum memory utilization percentage for the system',
      },
    },
    required: [],
  },
};

export async function handleUpdateSettings(
  args: {
    max_parallel_instances?: number;
    max_cpu_utilization_percent?: number;
    max_memory_utilization_percent?: number;
  },
  client: PySCFApiClient
) {
  try {
    // First get current settings
    const currentResponse = await client.getSettings();
    if (!currentResponse.success) {
      throw new Error('現在の設定取得に失敗しました');
    }

    const currentSettings = currentResponse.data.settings;

    // Merge with new settings, using safe defaults for undefined values
    const newSettings: SettingsUpdateRequest = {
      max_parallel_instances: args.max_parallel_instances ?? currentSettings.max_parallel_instances ?? 1,
      max_cpu_utilization_percent: args.max_cpu_utilization_percent ?? currentSettings.max_cpu_utilization_percent ?? 95,
      max_memory_utilization_percent: args.max_memory_utilization_percent ?? currentSettings.max_memory_utilization_percent ?? 95,
      system_total_cores: currentSettings.system_total_cores ?? 1,
      system_total_memory_mb: currentSettings.system_total_memory_mb ?? 1024,
    };

    const response = await client.updateSettings(newSettings);

    if (!response.success) {
      throw new Error(`設定の更新に失敗しました`);
    }

    const updatedSettings = response.data.settings;

    return {
      content: [
        {
          type: 'text',
          text: `✅ **設定を更新しました**

**更新された設定:**
- **最大並列インスタンス数:** ${updatedSettings.max_parallel_instances}
- **最大CPU使用率:** ${updatedSettings.max_cpu_utilization_percent}%
- **最大メモリ使用率:** ${updatedSettings.max_memory_utilization_percent}%

**実効制限:**
- **利用可能CPUコア:** ${Math.floor(updatedSettings.system_total_cores * updatedSettings.max_cpu_utilization_percent / 100)}コア
- **利用可能メモリ:** ${(updatedSettings.system_total_memory_mb * updatedSettings.max_memory_utilization_percent / 100 / 1024).toFixed(1)} GB

新しい設定は既に有効になっています。`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ 設定の更新でエラーが発生しました: ${errorMessage}

**可能な原因:**
- 無効な設定値
- サーバーエラー
- 権限不足

設定値が有効範囲内であることを確認してください。`,
        },
      ],
      isError: true,
    };
  }
}

export const getResourceStatusTool: Tool = {
  name: 'getResourceStatus',
  description: 'Get current system resource status including CPU and memory usage, and calculation resource allocation',
  inputSchema: {
    type: 'object',
    properties: {},
  },
};

export async function handleGetResourceStatus(_args: object, client: PySCFApiClient) {
  try {
    const response = await client.getResourceStatus();

    if (!response.success) {
      throw new Error(`リソース状況の取得に失敗しました`);
    }

    const data = response.data;
    const system = data.system_info;
    const constraints = data.resource_constraints;
    const allocated = data.allocated_resources;

    // Calculate percentages
    const cpuUsagePercent = system.cpu_usage_percent;
    const memoryUsagePercent = system.memory_usage_percent;
    const allocatedCpuPercent = (allocated.total_allocated_cpu_cores / system.total_cpu_cores * 100);
    const allocatedMemoryPercent = (allocated.total_allocated_memory_mb / system.total_memory_mb * 100);

    // Create progress bars
    const createProgressBar = (value: number, max: number = 100, width: number = 20) => {
      const percent = Math.min(value / max * 100, 100);
      const filled = Math.floor(percent / 100 * width);
      const empty = width - filled;
      return `[${'█'.repeat(filled)}${'░'.repeat(empty)}] ${percent.toFixed(1)}%`;
    };

    return {
      content: [
        {
          type: 'text',
          text: `📊 **システムリソース状況**

**現在のシステム使用状況:**
- **CPU使用率:** ${createProgressBar(cpuUsagePercent)} (${cpuUsagePercent.toFixed(1)}%)
- **メモリ使用率:** ${createProgressBar(memoryUsagePercent)} (${memoryUsagePercent.toFixed(1)}%)
- **利用可能メモリ:** ${(system.available_memory_mb / 1024).toFixed(1)} GB / ${(system.total_memory_mb / 1024).toFixed(1)} GB

**計算リソース制約:**
- **最大CPU利用率:** ${constraints.max_cpu_utilization_percent}%
- **最大メモリ利用率:** ${constraints.max_memory_utilization_percent}%
- **最大許可CPUコア:** ${constraints.max_allowed_cpu_cores}コア
- **最大許可メモリ:** ${(constraints.max_allowed_memory_mb / 1024).toFixed(1)} GB

**計算リソース割り当て:**
- **割り当て済みCPU:** ${createProgressBar(allocatedCpuPercent)} (${allocated.total_allocated_cpu_cores}/${system.total_cpu_cores}コア)
- **割り当て済みメモリ:** ${createProgressBar(allocatedMemoryPercent)} (${(allocated.total_allocated_memory_mb / 1024).toFixed(1)} GB)
- **利用可能CPUコア:** ${allocated.available_cpu_cores}コア
- **利用可能メモリ:** ${(allocated.available_memory_mb / 1024).toFixed(1)} GB
- **アクティブ計算数:** ${allocated.active_calculations_count}個

**更新時刻:** ${new Date(system.timestamp).toLocaleString()}

${allocated.active_calculations_count > 0 ? '現在実行中の計算があります。' : '現在実行中の計算はありません。'}`,
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
- レスポンスデータ: ${JSON.stringify(details.responseData, null, 2) || 'N/A'}

**考えられる原因:**
- サーバーが起動していない
- リソース管理モジュールの初期化エラー
- システム権限の問題
- ネットワーク接続の問題`;
    }
    
    return {
      content: [
        {
          type: 'text',
          text: `❌ リソース状況の取得でエラーが発生しました: ${errorMessage}${debugInfo}

**解決方法:**
- PySCF Native Appが正常に起動していることを確認してください
- \`testConnection\` ツールでサーバーとの接続を確認してください
- アプリケーションを再起動してください`,
        },
      ],
      isError: true,
    };
  }
}

export const testConnectionTool: Tool = {
  name: 'testConnection',
  description: 'Test connection to PySCF Native App server and get basic server information',
  inputSchema: {
    type: 'object',
    properties: {},
  },
};

export async function handleTestConnection(_args: object, client: PySCFApiClient) {
  try {
    const result = await client.testConnection();

    if (result.connected) {
      return {
        content: [
          {
            type: 'text',
            text: `✅ **PySCF Native Appサーバーに接続成功**

**サーバー情報:**
- **URL:** ${client.getBaseUrl()}
- **バージョン:** ${result.version || 'N/A'}
- **ステータス:** 正常稼働中

サーバーは正常に動作しており、全ての機能が利用可能です。`,
          },
        ],
      };
    } else {
      return {
        content: [
          {
            type: 'text',
            text: `❌ **PySCF Native Appサーバーに接続できません**

**エラー:** ${result.error}
**試行URL:** ${client.getBaseUrl()}

**解決方法:**
1. PySCF Native Appが起動していることを確認
2. \`npm run dev\` でアプリを起動
3. ポート番号が正しいことを確認 (通常5000-5100)
4. ファイアウォール設定を確認`,
          },
        ],
        isError: true,
      };
    }
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ 接続テストでエラーが発生しました: ${errorMessage}

サーバーの状態とネットワーク接続を確認してください。`,
        },
      ],
      isError: true,
    };
  }
}

export const diagnosticsServerTool: Tool = {
  name: 'diagnosticsServer',
  description: 'Perform comprehensive diagnostics of the PySCF Native App server including connectivity, endpoints, and dependencies',
  inputSchema: {
    type: 'object',
    properties: {
      detailed: {
        type: 'boolean',
        description: 'Include detailed endpoint testing and dependency information',
        default: false,
      },
    },
  },
};

export async function handleDiagnosticsServer(
  args: { detailed?: boolean },
  client: PySCFApiClient
) {
  const startTime = Date.now();
  const detailed = args.detailed || false;
  
  // Diagnostic results structure
  const results = {
    timestamp: new Date().toISOString(),
    server: { connected: false, baseUrl: client.getBaseUrl(), version: 'Unknown' },
    endpoints: { tested: 0, passed: 0, failed: 0, results: [] as any[] },
    dependencies: { available: [] as string[], missing: [] as string[] },
    recommendations: [] as string[]
  };

  let reportText = '🔍 **PySCF Native App サーバー診断レポート**\n\n';

  try {
    // 1. Basic connection test
    reportText += '**1. 基本接続テスト**\n';
    try {
      const connectionTest = await client.testConnection();
      results.server.connected = connectionTest.connected;
      results.server.version = connectionTest.version || 'Unknown';
      
      if (connectionTest.connected) {
        reportText += '✅ サーバー接続: 正常\n';
        reportText += `📍 サーバーURL: ${client.getBaseUrl()}\n`;
        reportText += `📋 サーバーバージョン: ${results.server.version}\n\n`;
      } else {
        reportText += '❌ サーバー接続: 失敗\n';
        reportText += `⚠️ エラー: ${connectionTest.error}\n\n`;
        results.recommendations.push('PySCF Native App が起動していることを確認してください');
        results.recommendations.push('npm run dev コマンドでアプリケーションを起動してください');
      }
    } catch (error) {
      reportText += '❌ サーバー接続: 失敗\n';
      reportText += `⚠️ エラー: ${error instanceof Error ? error.message : 'Unknown error'}\n\n`;
      results.recommendations.push('サーバーの起動状況を確認してください');
    }

    // 2. Critical endpoints testing (always performed)
    reportText += '**2. 重要エンドポイント診断**\n';
    
    // Test critical endpoints that were previously failing
    const criticalTests = [
      { name: 'getResourceStatus', test: () => client.getResourceStatus(), critical: true, description: '(以前404エラー)' },
      { name: 'getSupportedParameters', test: () => client.getSupportedParameters(), critical: false, description: '' },
    ];

    for (const endpoint of criticalTests) {
      results.endpoints.tested++;
      try {
        const result = await endpoint.test();
        results.endpoints.passed++;
        results.endpoints.results.push({ name: endpoint.name, status: 'PASS', error: null });
        reportText += `✅ ${endpoint.name}: 正常 ${endpoint.description}\n`;
        
        // Special handling for resource status to show CPU availability fix
        if (endpoint.name === 'getResourceStatus' && result.success) {
          // Type assertion for resource status response
          const resourceResult = result as any; // SystemResourceResponse type
          const data = resourceResult.data;
          if (data && data.allocated_resources && data.system_info) {
            const availableCpu = data.allocated_resources.available_cpu_cores;
            const totalCpu = data.system_info.total_cpu_cores;
            reportText += `   📊 CPU利用可能数: ${availableCpu}/${totalCpu}コア ${availableCpu > 0 ? '✅ 修正完了' : '❌ まだ問題'}\n`;
            reportText += `   🧠 メモリ利用可能: ${(data.allocated_resources.available_memory_mb / 1024).toFixed(1)} GB\n`;
          }
        }
      } catch (error) {
        results.endpoints.failed++;
        const errorMsg = error instanceof PySCFApiError ? 
          `${error.details.status} - ${error.details.statusText}` : 
          (error instanceof Error ? error.message : 'Unknown error');
        results.endpoints.results.push({ name: endpoint.name, status: 'FAIL', error: errorMsg });
        reportText += `❌ ${endpoint.name}: 失敗 ${endpoint.description} (${errorMsg})\n`;
        
        if (endpoint.name === 'getResourceStatus' && error instanceof PySCFApiError && error.details.status === 404) {
          results.recommendations.push('リソース管理モジュールの初期化に問題があります');
        }
      }
    }
    reportText += '\n';

    // 3. Detailed endpoint testing (only if detailed flag is set)
    if (results.server.connected && detailed) {
      reportText += '**3. 詳細エンドポイント診断**\n';
      
      const additionalTests = [
        { name: 'getSettings', test: () => client.getSettings() },
        { name: 'listCalculations', test: () => client.listCalculations() },
      ];

      for (const endpoint of additionalTests) {
        results.endpoints.tested++;
        try {
          await endpoint.test();
          results.endpoints.passed++;
          results.endpoints.results.push({ name: endpoint.name, status: 'PASS', error: null });
          reportText += `✅ ${endpoint.name}: 正常\n`;
        } catch (error) {
          results.endpoints.failed++;
          const errorMsg = error instanceof PySCFApiError ? 
            `${error.details.status} - ${error.details.statusText}` : 
            (error instanceof Error ? error.message : 'Unknown error');
          results.endpoints.results.push({ name: endpoint.name, status: 'FAIL', error: errorMsg });
          reportText += `❌ ${endpoint.name}: 失敗 (${errorMsg})\n`;
        }
      }
      reportText += '\n';
    }

    // 4. Parameter availability test
    if (results.server.connected) {
      reportText += '**4. パラメータ可用性テスト**\n';
      try {
        const params = await client.getSupportedParameters();
        if (params.success) {
          const data = params.data;
          reportText += `✅ 計算手法: ${data.calculation_methods.length}種類\n`;
          reportText += `✅ 基底関数: ${Object.values(data.basis_functions).flat().length}種類\n`;
          reportText += `✅ 交換相関汎関数: ${Object.values(data.exchange_correlation).flat().length}種類\n`;
          reportText += `✅ 溶媒: ${Object.values(data.solvents).flat().length}種類\n\n`;
          
          results.dependencies.available.push('計算パラメータ');
        }
      } catch (error) {
        reportText += '❌ パラメータ取得: 失敗\n';
        reportText += `⚠️ エラー: ${error instanceof Error ? error.message : 'Unknown error'}\n\n`;
        results.dependencies.missing.push('計算パラメータ');
        results.recommendations.push('量子化学ライブラリ（PySCF）の初期化を確認してください');
      }
    }

    // 5. StartCalculation test (only for detailed diagnostics)
    if (results.server.connected && detailed) {
      reportText += '**5. 計算開始機能テスト** (以前500エラー)\n';
      
      try {
        // Simple test molecule (hydrogen molecule)
        const testRequest = {
          xyz: `2
Hydrogen molecule test
H 0.0 0.0 0.0
H 0.74 0.0 0.0`,
          name: 'Diagnostic Test',
          calculation_method: 'DFT' as const,
          basis_function: 'STO-3G',
          exchange_correlation: 'B3LYP',
          charges: 0,
          spin: 0,
          solvent_method: 'none' as const,
          solvent: '-',
          cpu_cores: 1,
          memory_mb: 512,
          tddft_nstates: 5,
          tddft_method: 'TDDFT' as const,
          tddft_analyze_nto: false,
          // CASCI/CASSCF parameters (with defaults)
          ncas: 6,
          nelecas: 6,
          max_cycle_macro: 50,
          max_cycle_micro: 4,
          natorb: true,
          conv_tol: 0.000001,
          conv_tol_grad: 0.0001,
          optimize_geometry: false
        };
        
        const startResult = await client.startCalculation(testRequest);
        
        if (startResult.success) {
          const calc = startResult.data.calculation;
          reportText += `✅ 計算開始: 正常 (ID: ${calc.id})\n`;
          reportText += `   📋 ステータス: ${calc.status}\n`;
          
          // Try to immediately cancel the test calculation to avoid resource waste
          try {
            // Note: Cancel function not implemented in current client, so we skip this
            reportText += `   🔄 テスト計算のため、手動でキャンセルしてください\n`;
          } catch (e) {
            // Ignore cancel errors
          }
          
          results.dependencies.available.push('量子化学計算エンジン');
        } else {
          reportText += `❌ 計算開始: 失敗\n`;
          results.dependencies.missing.push('量子化学計算エンジン');
          results.recommendations.push('量子化学計算の初期化に問題があります');
        }
      } catch (error) {
        const errorMsg = error instanceof PySCFApiError ? 
          `${error.details.status} - ${error.details.statusText}` : 
          (error instanceof Error ? error.message : 'Unknown error');
        reportText += `❌ 計算開始: 失敗 (${errorMsg})\n`;
        
        if (error instanceof PySCFApiError && error.details.status === 500) {
          results.recommendations.push('サーバー内部エラー: 量子化学ライブラリまたはリソース管理に問題があります');
        }
        
        results.dependencies.missing.push('量子化学計算エンジン');
      }
      reportText += '\n';
    }

    // 6. Summary and recommendations
    const elapsedTime = Date.now() - startTime;
    reportText += '**📊 診断サマリー**\n';
    reportText += `⏱️ 診断時間: ${elapsedTime}ms\n`;
    reportText += `🔗 サーバー状態: ${results.server.connected ? '接続済み' : '未接続'}\n`;
    
    if (detailed) {
      reportText += `🎯 エンドポイントテスト: ${results.endpoints.passed}/${results.endpoints.tested} 成功\n`;
    }
    
    reportText += `📦 利用可能な依存関係: ${results.dependencies.available.length}\n`;
    reportText += `⚠️ 不足している依存関係: ${results.dependencies.missing.length}\n\n`;

    // 7. Recommendations
    if (results.recommendations.length > 0) {
      reportText += '**🔧 推奨事項**\n';
      results.recommendations.forEach((rec, i) => {
        reportText += `${i + 1}. ${rec}\n`;
      });
      reportText += '\n';
    }

    // 8. Overall health status
    const overallHealth = results.server.connected && results.endpoints.failed === 0;
    reportText += `**🏥 総合ヘルス状態: ${overallHealth ? '健全' : '要注意'}**\n`;
    
    if (overallHealth) {
      reportText += '✅ システムは正常に動作しています。';
    } else {
      reportText += '⚠️ 一部機能に問題があります。上記の推奨事項をご確認ください。';
    }

    return {
      content: [
        {
          type: 'text',
          text: reportText,
        },
      ],
    };

  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ 診断実行でエラーが発生しました: ${errorMessage}

**問題:**
診断プロセス自体でエラーが発生しました。

**解決方法:**
1. ネットワーク接続を確認してください
2. PySCF Native App の起動状況を確認してください
3. しばらく待ってから再試行してください`,
        },
      ],
      isError: true,
    };
  }
}