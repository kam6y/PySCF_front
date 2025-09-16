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
    results.push('🔍 **Diagnostic Results Report**\n');
    
    const connectionResult = await client.testConnection();
    const connectionTime = Date.now() - startTime;
    
    if (connectionResult.connected) {
      results.push(`✅ **Connection Test**: Success (${connectionTime}ms)`);
      results.push(`   - Server URL: ${client.getBaseUrl()}`);
      results.push(`   - Version: ${connectionResult.version || 'N/A'}`);
    } else {
      results.push(`❌ **Connection Test**: Failed`);
      results.push(`   - Error: ${connectionResult.error}`);
      
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
        results.push(`   - Health Check: Normal`);
        results.push(`   - Service: ${healthData.service || 'N/A'}`);
        results.push(`   - Status: ${healthData.status || 'N/A'}`);
      } catch (error) {
        results.push(`   - Health Check: Error (${error})`);
      }
    }

    // 3. Endpoint Availability
    results.push('\n📡 **Endpoint Availability**');
    const endpoints = [
      { name: 'Supported Parameters', method: 'getSupportedParameters' },
      { name: 'Get Settings', method: 'getSettings' },
      { name: 'Resource Status', method: 'getResourceStatus' },
      { name: 'Calculation List', method: 'listCalculations' },
    ];

    for (const endpoint of endpoints) {
      try {
        const testStart = Date.now();
        await (client as any)[endpoint.method]();
        const testTime = Date.now() - testStart;
        results.push(`✅ ${endpoint.name}: Available (${testTime}ms)`);
      } catch (error) {
        results.push(`❌ ${endpoint.name}: Error (${error})`);
      }
    }

    // 4. Performance Metrics
    if (args.include_performance) {
      results.push('\n📊 **Performance Statistics**');
      const stats = MCPInspector.getErrorStats();
      results.push(`- Total Requests: ${stats.total}`);
      results.push(`- Errors: ${stats.errors}`);
      results.push(`- Success Rate: ${stats.successRate}%`);
      
      const recentLogs = MCPInspector.getLogs(10);
      if (recentLogs.length > 0) {
        const avgResponseTime = recentLogs
          .filter(log => log.responseTime)
          .reduce((sum, log) => sum + (log.responseTime || 0), 0) / recentLogs.length;
        results.push(`- Average Response Time: ${avgResponseTime.toFixed(1)}ms`);
      }
    }

    // 5. Recent Errors
    const recentErrors = MCPInspector.getLogs(5).filter(log => !log.success);
    if (recentErrors.length > 0) {
      results.push('\n⚠️ **Recent Errors**');
      recentErrors.forEach((log, index) => {
        results.push(`${index + 1}. [${log.timestamp.toLocaleTimeString()}] ${log.tool}: ${log.error}`);
      });
    }

    const totalTime = Date.now() - startTime;
    results.push(`\nDiagnostics completed (total execution time: ${totalTime}ms)`);

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
          text: `❌ Error occurred during diagnostics: ${errorMessage}`,
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

  results.push(`🧪 **API Endpoint Test** (${endpoint})\n`);

  const testEndpoint = async (name: string, testFunc: () => Promise<any>) => {
    try {
      const testStart = Date.now();
      const response = await testFunc();
      const responseTime = Date.now() - testStart;
      
      results.push(`✅ **${name}**: Success (${responseTime}ms)`);
      
      if (args.include_response_data && response) {
        const dataStructure = typeof response === 'object' 
          ? Object.keys(response).join(', ')
          : typeof response;
        results.push(`   - Response structure: ${dataStructure}`);
      }
      
      return true;
    } catch (error) {
      const errorMessage = error instanceof Error ? error.message : 'Unknown error';
      results.push(`❌ **${name}**: Error - ${errorMessage}`);
      return false;
    }
  };

  try {
    if (endpoint === 'all' || endpoint === 'system') {
      await testEndpoint('Health Check', () => client.healthCheck());
      await testEndpoint('Get Settings', () => client.getSettings());
      await testEndpoint('Resource Status', () => client.getResourceStatus());
      await testEndpoint('Supported Parameters', () => client.getSupportedParameters());
    }

    if (endpoint === 'all' || endpoint === 'quantum') {
      await testEndpoint('Calculation List', () => client.listCalculations());
    }

    if (endpoint === 'all' || endpoint === 'pubchem') {
      await testEndpoint('XYZ Validation', () =>
        client.validateXYZ({ xyz: '1\nTest\nH 0.0 0.0 0.0' })
      );
    }

    const totalTime = Date.now() - startTime;
    results.push(`\nTest completed (total execution time: ${totalTime}ms)`);

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
          text: `❌ Error occurred during API testing: ${errorMessage}`,
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
            text: `📝 **Debug Logs**\n\n${args.errors_only ? 'No error logs available.' : 'No logs available.'}\n\nOperation history will be displayed as it becomes available.`,
          },
        ],
      };
    }

    const results: string[] = [];
    results.push(`📝 **Debug Logs** (Latest ${logs.length} entries)\n`);

    logs.forEach((log, index) => {
      const status = log.success ? '✅' : '❌';
      const time = log.timestamp.toLocaleTimeString();
      const responseTime = log.responseTime ? ` (${log.responseTime}ms)` : '';
      
      results.push(`${index + 1}. ${status} [${time}] **${log.tool}**${responseTime}`);
      
      if (log.args && Object.keys(log.args).length > 0) {
        const argStr = Object.entries(log.args)
          .map(([key, value]) => `${key}: ${JSON.stringify(value)}`)
          .join(', ');
        results.push(`   - Arguments: {${argStr}}`);
      }
      
      if (log.error) {
        results.push(`   - Error: ${log.error}`);
      }
      
      results.push('');
    });

    // Error pattern analysis
    const errorCount = logs.filter(log => !log.success).length;
    if (errorCount > 0) {
      results.push('🔍 **Error Pattern Analysis**');
      results.push(`- Total Errors: ${errorCount} times`);
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
          text: `❌ Error occurred while retrieving logs: ${errorMessage}`,
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
    results.push('🔧 **MCP Configuration Validation Report**\n');

    // 1. Connection Configuration
    results.push('**Connection Settings**');
    const baseUrl = client.getBaseUrl();
    results.push(`- Server URL: ${baseUrl}`);
    
    if (baseUrl.includes('127.0.0.1') || baseUrl.includes('localhost')) {
      results.push('  ✅ Appropriate for local development environment');
    } else {
      results.push('  ⚠️ For remote servers, please verify security settings');
    }

    // 2. Server Settings Validation
    try {
      const settingsResponse = await client.getSettings();
      if (settingsResponse.success) {
        const settings = settingsResponse.data.settings;
        results.push('\n**Server Settings Validation**');
        
        // Validate parallel instances
        const maxInstances = settings.max_parallel_instances || 1;
        if (maxInstances >= 1 && maxInstances <= 8) {
          results.push(`✅ Parallel instances: ${maxInstances} (within recommended range)`);
        } else if (maxInstances > 8) {
          results.push(`⚠️ Parallel instances: ${maxInstances} (may be too many)`);
        } else {
          results.push(`❌ Parallel instances: ${maxInstances} (invalid value)`);
        }
        
        // Validate resource limits
        const cpuLimit = settings.max_cpu_utilization_percent || 95;
        const memLimit = settings.max_memory_utilization_percent || 95;
        
        if (cpuLimit >= 50 && cpuLimit <= 95) {
          results.push(`✅ CPU usage limit: ${cpuLimit}% (appropriate)`);
        } else {
          results.push(`⚠️ CPU usage limit: ${cpuLimit}% (50-95% recommended)`);
        }
        
        if (memLimit >= 50 && memLimit <= 95) {
          results.push(`✅ Memory usage limit: ${memLimit}% (appropriate)`);
        } else {
          results.push(`⚠️ Memory usage limit: ${memLimit}% (50-95% recommended)`);
        }
      }
    } catch (error) {
      results.push('\n❌ **Failed to retrieve server settings**');
      results.push(`Error: ${error}`);
    }

    // 3. Performance Recommendations
    if (args.check_performance) {
      results.push('\n**Performance Optimization Suggestions**');
      
      const stats = MCPInspector.getErrorStats();
      if (parseInt(stats.successRate) < 90) {
        results.push('⚠️ Success rate is below 90%');
        results.push('  - Check network connection');
        results.push('  - Check server resources');
      } else {
        results.push('✅ API success rate is good');
      }
      
      // Response time analysis
      const recentLogs = MCPInspector.getLogs(20);
      const responseTimes = recentLogs
        .filter(log => log.responseTime && log.responseTime > 0)
        .map(log => log.responseTime!);
      
      if (responseTimes.length > 0) {
        const avgTime = responseTimes.reduce((a, b) => a + b, 0) / responseTimes.length;
        if (avgTime > 5000) {
          results.push('⚠️ Average response time exceeds 5 seconds');
          results.push('  - Check timeout settings');
          results.push('  - Check server load');
        } else {
          results.push('✅ Response time is good');
        }
      }
    }

    // 4. Recommendations
    results.push('\n**Recommendations**');
    results.push('1. Run `diagnostics` regularly to check server status');
    results.push('2. If errors occur frequently, analyze patterns with `getDebugLogs`');
    results.push('3. Check resources with `getResourceStatus` before heavy calculations');
    results.push('4. Verify available parameters with `getSupportedParameters`');

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
          text: `❌ Error occurred during configuration validation: ${errorMessage}`,
        },
      ],
      isError: true,
    };
  }
}

// Export the inspector class for potential use in other tools
export { MCPInspector };