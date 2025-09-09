import { Tool } from '@modelcontextprotocol/sdk/types.js';
import { PySCFApiClient } from '../client.js';

// Debug log storage
interface DebugLog {
  timestamp: Date;
  tool: string;
  args: any;
  success: boolean;
  error?: string;
  responseTime?: number;
}

class MCPInspector {
  private static logs: DebugLog[] = [];
  private static readonly MAX_LOGS = 100;

  static addLog(log: DebugLog) {
    this.logs.unshift(log);
    if (this.logs.length > this.MAX_LOGS) {
      this.logs = this.logs.slice(0, this.MAX_LOGS);
    }
  }

  static getLogs(limit?: number): DebugLog[] {
    return limit ? this.logs.slice(0, limit) : this.logs;
  }

  static getErrorStats() {
    const total = this.logs.length;
    const errors = this.logs.filter(log => !log.success).length;
    return { total, errors, successRate: total > 0 ? ((total - errors) / total * 100).toFixed(1) : '0' };
  }
}

export const diagnosticsTool: Tool = {
  name: 'diagnostics',
  description: 'Run comprehensive diagnostics on the PySCF MCP server and connection',
  inputSchema: {
    type: 'object',
    properties: {
      include_performance: {
        type: 'boolean',
        description: 'Include performance metrics in the report',
        default: true,
      },
    },
  },
};

export async function handleDiagnostics(
  args: { include_performance?: boolean },
  client: PySCFApiClient
) {
  const startTime = Date.now();
  const results: string[] = [];

  try {
    // 1. Connection Test
    results.push('🔍 **診断結果レポート**\n');
    
    const connectionResult = await client.testConnection();
    const connectionTime = Date.now() - startTime;
    
    if (connectionResult.connected) {
      results.push(`✅ **接続テスト**: 成功 (${connectionTime}ms)`);
      results.push(`   - サーバーURL: ${client.getBaseUrl()}`);
      results.push(`   - バージョン: ${connectionResult.version || 'N/A'}`);
    } else {
      results.push(`❌ **接続テスト**: 失敗`);
      results.push(`   - エラー: ${connectionResult.error}`);
      
      MCPInspector.addLog({
        timestamp: new Date(),
        tool: 'diagnostics',
        args,
        success: false,
        error: connectionResult.error || 'Connection failed',
        responseTime: connectionTime,
      });
      
      return {
        content: [{ type: 'text', text: results.join('\n') }],
        isError: true,
      };
    }

    // 2. Health Check Details
    if (connectionResult.connected) {
      try {
        const healthData = await client.healthCheck();
        results.push(`   - ヘルスチェック: 正常`);
        results.push(`   - サービス: ${healthData.service || 'N/A'}`);
        results.push(`   - ステータス: ${healthData.status || 'N/A'}`);
      } catch (error) {
        results.push(`   - ヘルスチェック: エラー (${error})`);
      }
    }

    // 3. Endpoint Availability
    results.push('\n📡 **エンドポイント可用性**');
    const endpoints = [
      { name: 'サポートパラメータ', method: 'getSupportedParameters' },
      { name: '設定取得', method: 'getSettings' },
      { name: 'リソース状況', method: 'getResourceStatus' },
      { name: '計算リスト', method: 'listCalculations' },
    ];

    for (const endpoint of endpoints) {
      try {
        const testStart = Date.now();
        await (client as any)[endpoint.method]();
        const testTime = Date.now() - testStart;
        results.push(`✅ ${endpoint.name}: 利用可能 (${testTime}ms)`);
      } catch (error) {
        results.push(`❌ ${endpoint.name}: エラー (${error})`);
      }
    }

    // 4. Performance Metrics
    if (args.include_performance) {
      results.push('\n📊 **パフォーマンス統計**');
      const stats = MCPInspector.getErrorStats();
      results.push(`- 総リクエスト数: ${stats.total}`);
      results.push(`- エラー数: ${stats.errors}`);
      results.push(`- 成功率: ${stats.successRate}%`);
      
      const recentLogs = MCPInspector.getLogs(10);
      if (recentLogs.length > 0) {
        const avgResponseTime = recentLogs
          .filter(log => log.responseTime)
          .reduce((sum, log) => sum + (log.responseTime || 0), 0) / recentLogs.length;
        results.push(`- 平均応答時間: ${avgResponseTime.toFixed(1)}ms`);
      }
    }

    // 5. Recent Errors
    const recentErrors = MCPInspector.getLogs(5).filter(log => !log.success);
    if (recentErrors.length > 0) {
      results.push('\n⚠️ **最近のエラー**');
      recentErrors.forEach((log, index) => {
        results.push(`${index + 1}. [${log.timestamp.toLocaleTimeString()}] ${log.tool}: ${log.error}`);
      });
    }

    const totalTime = Date.now() - startTime;
    results.push(`\n診断完了 (総実行時間: ${totalTime}ms)`);

    MCPInspector.addLog({
      timestamp: new Date(),
      tool: 'diagnostics',
      args,
      success: true,
      responseTime: totalTime,
    });

    return {
      content: [{ type: 'text', text: results.join('\n') }],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error';
    
    MCPInspector.addLog({
      timestamp: new Date(),
      tool: 'diagnostics',
      args,
      success: false,
      error: errorMessage,
      responseTime: Date.now() - startTime,
    });

    return {
      content: [
        {
          type: 'text',
          text: `❌ 診断中にエラーが発生しました: ${errorMessage}`,
        },
      ],
      isError: true,
    };
  }
}

export const testApiTool: Tool = {
  name: 'testApi',
  description: 'Test specific API endpoints with response time and validation',
  inputSchema: {
    type: 'object',
    properties: {
      endpoint: {
        type: 'string',
        enum: ['all', 'pubchem', 'quantum', 'system'],
        description: 'Which API endpoints to test',
        default: 'all',
      },
      include_response_data: {
        type: 'boolean',
        description: 'Include response data structure in the output',
        default: false,
      },
    },
  },
};

export async function handleTestApi(
  args: { endpoint?: string; include_response_data?: boolean },
  client: PySCFApiClient
) {
  const startTime = Date.now();
  const results: string[] = [];
  const endpoint = args.endpoint || 'all';

  results.push(`🧪 **API エンドポイントテスト** (${endpoint})\n`);

  const testEndpoint = async (name: string, testFunc: () => Promise<any>) => {
    try {
      const testStart = Date.now();
      const response = await testFunc();
      const responseTime = Date.now() - testStart;
      
      results.push(`✅ **${name}**: 成功 (${responseTime}ms)`);
      
      if (args.include_response_data && response) {
        const dataStructure = typeof response === 'object' 
          ? Object.keys(response).join(', ')
          : typeof response;
        results.push(`   - レスポンス構造: ${dataStructure}`);
      }
      
      return true;
    } catch (error) {
      const errorMessage = error instanceof Error ? error.message : 'Unknown error';
      results.push(`❌ **${name}**: エラー - ${errorMessage}`);
      return false;
    }
  };

  try {
    if (endpoint === 'all' || endpoint === 'system') {
      await testEndpoint('ヘルスチェック', () => client.healthCheck());
      await testEndpoint('設定取得', () => client.getSettings());
      await testEndpoint('リソース状況', () => client.getResourceStatus());
      await testEndpoint('サポートパラメータ', () => client.getSupportedParameters());
    }

    if (endpoint === 'all' || endpoint === 'quantum') {
      await testEndpoint('計算リスト', () => client.listCalculations());
    }

    if (endpoint === 'all' || endpoint === 'pubchem') {
      await testEndpoint('XYZ検証', () => 
        client.validateXYZ({ xyz: '1\nTest\nH 0.0 0.0 0.0' })
      );
    }

    const totalTime = Date.now() - startTime;
    results.push(`\nテスト完了 (総実行時間: ${totalTime}ms)`);

    MCPInspector.addLog({
      timestamp: new Date(),
      tool: 'testApi',
      args,
      success: true,
      responseTime: totalTime,
    });

    return {
      content: [{ type: 'text', text: results.join('\n') }],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error';
    
    MCPInspector.addLog({
      timestamp: new Date(),
      tool: 'testApi',
      args,
      success: false,
      error: errorMessage,
      responseTime: Date.now() - startTime,
    });

    return {
      content: [
        {
          type: 'text',
          text: `❌ APIテスト中にエラーが発生しました: ${errorMessage}`,
        },
      ],
      isError: true,
    };
  }
}

export const getDebugLogsTool: Tool = {
  name: 'getDebugLogs',
  description: 'Get recent debug logs and error patterns from MCP operations',
  inputSchema: {
    type: 'object',
    properties: {
      limit: {
        type: 'integer',
        minimum: 1,
        maximum: 50,
        description: 'Number of recent logs to retrieve',
        default: 20,
      },
      errors_only: {
        type: 'boolean',
        description: 'Show only error logs',
        default: false,
      },
    },
  },
};

export async function handleGetDebugLogs(
  args: { limit?: number; errors_only?: boolean }
) {
  try {
    const limit = args.limit || 20;
    let logs = MCPInspector.getLogs(limit);
    
    if (args.errors_only) {
      logs = logs.filter(log => !log.success);
    }

    if (logs.length === 0) {
      return {
        content: [
          {
            type: 'text',
            text: `📝 **デバッグログ**\n\n${args.errors_only ? 'エラーログがありません。' : 'ログがありません。'}\n\n使用可能になり次第、操作履歴が表示されます。`,
          },
        ],
      };
    }

    const results: string[] = [];
    results.push(`📝 **デバッグログ** (最新 ${logs.length} 件)\n`);

    logs.forEach((log, index) => {
      const status = log.success ? '✅' : '❌';
      const time = log.timestamp.toLocaleTimeString();
      const responseTime = log.responseTime ? ` (${log.responseTime}ms)` : '';
      
      results.push(`${index + 1}. ${status} [${time}] **${log.tool}**${responseTime}`);
      
      if (log.args && Object.keys(log.args).length > 0) {
        const argStr = Object.entries(log.args)
          .map(([key, value]) => `${key}: ${JSON.stringify(value)}`)
          .join(', ');
        results.push(`   - 引数: {${argStr}}`);
      }
      
      if (log.error) {
        results.push(`   - エラー: ${log.error}`);
      }
      
      results.push('');
    });

    // Error pattern analysis
    const errorCount = logs.filter(log => !log.success).length;
    if (errorCount > 0) {
      results.push('🔍 **エラーパターン分析**');
      results.push(`- 総エラー数: ${errorCount}回`);
    }

    return {
      content: [{ type: 'text', text: results.join('\n') }],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error';
    return {
      content: [
        {
          type: 'text',
          text: `❌ ログ取得でエラーが発生しました: ${errorMessage}`,
        },
      ],
      isError: true,
    };
  }
}

export const validateConfigTool: Tool = {
  name: 'validateConfig',
  description: 'Validate MCP configuration and provide optimization recommendations',
  inputSchema: {
    type: 'object',
    properties: {
      check_performance: {
        type: 'boolean',
        description: 'Include performance optimization suggestions',
        default: true,
      },
    },
  },
};

export async function handleValidateConfig(
  args: { check_performance?: boolean },
  client: PySCFApiClient
) {
  try {
    const results: string[] = [];
    results.push('🔧 **MCP設定検証レポート**\n');

    // 1. Connection Configuration
    results.push('**接続設定**');
    const baseUrl = client.getBaseUrl();
    results.push(`- サーバーURL: ${baseUrl}`);
    
    if (baseUrl.includes('127.0.0.1') || baseUrl.includes('localhost')) {
      results.push('  ✅ ローカル開発環境として適切');
    } else {
      results.push('  ⚠️ リモートサーバーの場合、セキュリティ設定を確認してください');
    }

    // 2. Server Settings Validation
    try {
      const settingsResponse = await client.getSettings();
      if (settingsResponse.success) {
        const settings = settingsResponse.data.settings;
        results.push('\n**サーバー設定検証**');
        
        // Validate parallel instances
        const maxInstances = settings.max_parallel_instances || 1;
        if (maxInstances >= 1 && maxInstances <= 8) {
          results.push(`✅ 並列インスタンス数: ${maxInstances} (推奨範囲内)`);
        } else if (maxInstances > 8) {
          results.push(`⚠️ 並列インスタンス数: ${maxInstances} (多すぎる可能性があります)`);
        } else {
          results.push(`❌ 並列インスタンス数: ${maxInstances} (無効な値)`);
        }
        
        // Validate resource limits
        const cpuLimit = settings.max_cpu_utilization_percent || 95;
        const memLimit = settings.max_memory_utilization_percent || 95;
        
        if (cpuLimit >= 50 && cpuLimit <= 95) {
          results.push(`✅ CPU使用率制限: ${cpuLimit}% (適切)`);
        } else {
          results.push(`⚠️ CPU使用率制限: ${cpuLimit}% (50-95%を推奨)`);
        }
        
        if (memLimit >= 50 && memLimit <= 95) {
          results.push(`✅ メモリ使用率制限: ${memLimit}% (適切)`);
        } else {
          results.push(`⚠️ メモリ使用率制限: ${memLimit}% (50-95%を推奨)`);
        }
      }
    } catch (error) {
      results.push('\n❌ **サーバー設定の取得に失敗**');
      results.push(`エラー: ${error}`);
    }

    // 3. Performance Recommendations
    if (args.check_performance) {
      results.push('\n**パフォーマンス最適化提案**');
      
      const stats = MCPInspector.getErrorStats();
      if (parseInt(stats.successRate) < 90) {
        results.push('⚠️ 成功率が90%を下回っています');
        results.push('  - ネットワーク接続を確認してください');
        results.push('  - サーバーリソースを確認してください');
      } else {
        results.push('✅ API成功率は良好です');
      }
      
      // Response time analysis
      const recentLogs = MCPInspector.getLogs(20);
      const responseTimes = recentLogs
        .filter(log => log.responseTime && log.responseTime > 0)
        .map(log => log.responseTime!);
      
      if (responseTimes.length > 0) {
        const avgTime = responseTimes.reduce((a, b) => a + b, 0) / responseTimes.length;
        if (avgTime > 5000) {
          results.push('⚠️ 平均応答時間が5秒を超えています');
          results.push('  - タイムアウト設定を確認してください');
          results.push('  - サーバーの負荷を確認してください');
        } else {
          results.push('✅ 応答時間は良好です');
        }
      }
    }

    // 4. Recommendations
    results.push('\n**推奨事項**');
    results.push('1. 定期的に`diagnostics`を実行してサーバー状況を確認');
    results.push('2. エラーが頻発する場合は`getDebugLogs`でパターンを分析');
    results.push('3. 重い計算の前に`getResourceStatus`でリソースを確認');
    results.push('4. `getSupportedParameters`で利用可能なパラメータを確認');

    MCPInspector.addLog({
      timestamp: new Date(),
      tool: 'validateConfig',
      args,
      success: true,
    });

    return {
      content: [{ type: 'text', text: results.join('\n') }],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error';
    
    MCPInspector.addLog({
      timestamp: new Date(),
      tool: 'validateConfig',
      args,
      success: false,
      error: errorMessage,
    });

    return {
      content: [
        {
          type: 'text',
          text: `❌ 設定検証中にエラーが発生しました: ${errorMessage}`,
        },
      ],
      isError: true,
    };
  }
}

// Export the inspector class for potential use in other tools
export { MCPInspector };