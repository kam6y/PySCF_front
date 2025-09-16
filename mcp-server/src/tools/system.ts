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
      throw new Error(`Failed to get supported parameters`);
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
          text: `⚙️ **Supported Calculation Parameters**

**Calculation Methods:**
${data.calculation_methods.map(method => `- ${method}`).join('\n')}

**Basis Functions:**
${basisText}

**Exchange-Correlation Functionals:**
${xcText}

**Solvent Effect Methods:**
${data.solvent_methods.map(method => `- ${method}`).join('\n')}

**Solvents:**
${solventText}

**TDDFT Methods:**
${data.tddft_methods.map(method => `- ${method}`).join('\n')}

These parameters can be used with \`startCalculation\`.`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ Error occurred while getting supported parameters: ${errorMessage}

Please ensure the server is running.`,
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
      throw new Error(`Failed to get settings`);
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
          text: `⚙️ **Application Settings**

**Parallel Processing Settings:**
- **Max Parallel Instances:** ${safeSettings.max_parallel_instances}
- **Max CPU Utilization:** ${safeSettings.max_cpu_utilization_percent}%
- **Max Memory Utilization:** ${safeSettings.max_memory_utilization_percent}%

**System Information:**
- **Total CPU Cores:** ${safeSettings.system_total_cores}
- **Total Memory:** ${(safeSettings.system_total_memory_mb / 1024).toFixed(1)} GB (${safeSettings.system_total_memory_mb} MB)

**Effective Limits:**
- **Available CPU Cores:** ${Math.floor(safeSettings.system_total_cores * safeSettings.max_cpu_utilization_percent / 100)} cores
- **Available Memory:** ${(safeSettings.system_total_memory_mb * safeSettings.max_memory_utilization_percent / 100 / 1024).toFixed(1)} GB

Use \`updateSettings\` to change these settings.`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ Error occurred while getting settings: ${errorMessage}`,
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
      throw new Error('Failed to get current settings');
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
      throw new Error(`Failed to update settings`);
    }

    const updatedSettings = response.data.settings;

    return {
      content: [
        {
          type: 'text',
          text: `✅ **Settings Updated Successfully**

**Updated Settings:**
- **Max Parallel Instances:** ${updatedSettings.max_parallel_instances}
- **Max CPU Utilization:** ${updatedSettings.max_cpu_utilization_percent}%
- **Max Memory Utilization:** ${updatedSettings.max_memory_utilization_percent}%

**Effective Limits:**
- **Available CPU Cores:** ${Math.floor(updatedSettings.system_total_cores * updatedSettings.max_cpu_utilization_percent / 100)} cores
- **Available Memory:** ${(updatedSettings.system_total_memory_mb * updatedSettings.max_memory_utilization_percent / 100 / 1024).toFixed(1)} GB

The new settings are now active.`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ Error occurred while updating settings: ${errorMessage}

**Possible Causes:**
- Invalid setting values
- Server error
- Insufficient permissions

Please ensure the setting values are within valid ranges.`,
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
      throw new Error(`Failed to get resource status`);
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
          text: `📊 **System Resource Status**

**Current System Usage:**
- **CPU Usage:** ${createProgressBar(cpuUsagePercent)} (${cpuUsagePercent.toFixed(1)}%)
- **Memory Usage:** ${createProgressBar(memoryUsagePercent)} (${memoryUsagePercent.toFixed(1)}%)
- **Available Memory:** ${(system.available_memory_mb / 1024).toFixed(1)} GB / ${(system.total_memory_mb / 1024).toFixed(1)} GB

**Calculation Resource Constraints:**
- **Max CPU Utilization:** ${constraints.max_cpu_utilization_percent}%
- **Max Memory Utilization:** ${constraints.max_memory_utilization_percent}%
- **Max Allowed CPU Cores:** ${constraints.max_allowed_cpu_cores} cores
- **Max Allowed Memory:** ${(constraints.max_allowed_memory_mb / 1024).toFixed(1)} GB

**Calculation Resource Allocation:**
- **Allocated CPU:** ${createProgressBar(allocatedCpuPercent)} (${allocated.total_allocated_cpu_cores}/${system.total_cpu_cores} cores)
- **Allocated Memory:** ${createProgressBar(allocatedMemoryPercent)} (${(allocated.total_allocated_memory_mb / 1024).toFixed(1)} GB)
- **Available CPU Cores:** ${allocated.available_cpu_cores} cores
- **Available Memory:** ${(allocated.available_memory_mb / 1024).toFixed(1)} GB
- **Active Calculations:** ${allocated.active_calculations_count}

**Last Updated:** ${new Date(system.timestamp).toLocaleString()}

${allocated.active_calculations_count > 0 ? 'There are currently running calculations.' : 'No calculations are currently running.'}`,
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
- Response Data: ${JSON.stringify(details.responseData, null, 2) || 'N/A'}

**Possible Causes:**
- Server is not running
- Resource management module initialization error
- System permission issues
- Network connection problems`;
    }
    
    return {
      content: [
        {
          type: 'text',
          text: `❌ Error occurred while getting resource status: ${errorMessage}${debugInfo}

**Solutions:**
- Ensure PySCF Native App is running properly
- Use \`testConnection\` tool to verify server connection
- Restart the application`,
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
            text: `✅ **Successfully connected to PySCF Native App server**

**Server Information:**
- **URL:** ${client.getBaseUrl()}
- **Version:** ${result.version || 'N/A'}
- **Status:** Running normally

The server is operating properly and all features are available.`,
          },
        ],
      };
    } else {
      return {
        content: [
          {
            type: 'text',
            text: `❌ **Cannot connect to PySCF Native App server**

**Error:** ${result.error}
**Attempted URL:** ${client.getBaseUrl()}

**Solutions:**
1. Ensure PySCF Native App is running
2. Start the app with \`npm run dev\`
3. Verify the port number is correct (usually 5000-5100)
4. Check firewall settings`,
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
          text: `❌ Error occurred during connection test: ${errorMessage}

Please check server status and network connection.`,
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

  let reportText = '🔍 **PySCF Native App Server Diagnostic Report**\n\n';

  try {
    // 1. Basic connection test
    reportText += '**1. Basic Connection Test**\n';
    try {
      const connectionTest = await client.testConnection();
      results.server.connected = connectionTest.connected;
      results.server.version = connectionTest.version || 'Unknown';
      
      if (connectionTest.connected) {
        reportText += '✅ Server connection: Normal\n';
        reportText += `📍 Server URL: ${client.getBaseUrl()}\n`;
        reportText += `📋 Server version: ${results.server.version}\n\n`;
      } else {
        reportText += '❌ Server connection: Failed\n';
        reportText += `⚠️ Error: ${connectionTest.error}\n\n`;
        results.recommendations.push('Ensure PySCF Native App is running');
        results.recommendations.push('Start the application with npm run dev command');
      }
    } catch (error) {
      reportText += '❌ Server connection: Failed\n';
      reportText += `⚠️ Error: ${error instanceof Error ? error.message : 'Unknown error'}\n\n`;
      results.recommendations.push('Check server startup status');
    }

    // 2. Critical endpoints testing (always performed)
    reportText += '**2. Critical Endpoint Diagnostics**\n';
    
    // Test critical endpoints that were previously failing
    const criticalTests = [
      { name: 'getResourceStatus', test: () => client.getResourceStatus(), critical: true, description: '(previously 404 error)' },
      { name: 'getSupportedParameters', test: () => client.getSupportedParameters(), critical: false, description: '' },
    ];

    for (const endpoint of criticalTests) {
      results.endpoints.tested++;
      try {
        const result = await endpoint.test();
        results.endpoints.passed++;
        results.endpoints.results.push({ name: endpoint.name, status: 'PASS', error: null });
        reportText += `✅ ${endpoint.name}: Normal ${endpoint.description}\n`;
        
        // Special handling for resource status to show CPU availability fix
        if (endpoint.name === 'getResourceStatus' && result.success) {
          // Type assertion for resource status response
          const resourceResult = result as any; // SystemResourceResponse type
          const data = resourceResult.data;
          if (data && data.allocated_resources && data.system_info) {
            const availableCpu = data.allocated_resources.available_cpu_cores;
            const totalCpu = data.system_info.total_cpu_cores;
            reportText += `   📊 Available CPU cores: ${availableCpu}/${totalCpu} cores ${availableCpu > 0 ? '✅ Fixed' : '❌ Still problematic'}\n`;
            reportText += `   🧠 Available memory: ${(data.allocated_resources.available_memory_mb / 1024).toFixed(1)} GB\n`;
          }
        }
      } catch (error) {
        results.endpoints.failed++;
        const errorMsg = error instanceof PySCFApiError ? 
          `${error.details.status} - ${error.details.statusText}` : 
          (error instanceof Error ? error.message : 'Unknown error');
        results.endpoints.results.push({ name: endpoint.name, status: 'FAIL', error: errorMsg });
        reportText += `❌ ${endpoint.name}: Failed ${endpoint.description} (${errorMsg})\n`;
        
        if (endpoint.name === 'getResourceStatus' && error instanceof PySCFApiError && error.details.status === 404) {
          results.recommendations.push('Resource management module initialization issue');
        }
      }
    }
    reportText += '\n';

    // 3. Detailed endpoint testing (only if detailed flag is set)
    if (results.server.connected && detailed) {
      reportText += '**3. Detailed Endpoint Diagnostics**\n';
      
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
          reportText += `✅ ${endpoint.name}: Normal\n`;
        } catch (error) {
          results.endpoints.failed++;
          const errorMsg = error instanceof PySCFApiError ? 
            `${error.details.status} - ${error.details.statusText}` : 
            (error instanceof Error ? error.message : 'Unknown error');
          results.endpoints.results.push({ name: endpoint.name, status: 'FAIL', error: errorMsg });
          reportText += `❌ ${endpoint.name}: Failed (${errorMsg})\n`;
        }
      }
      reportText += '\n';
    }

    // 4. Parameter availability test
    if (results.server.connected) {
      reportText += '**4. Parameter Availability Test**\n';
      try {
        const params = await client.getSupportedParameters();
        if (params.success) {
          const data = params.data;
          reportText += `✅ Calculation methods: ${data.calculation_methods.length} types\n`;
          reportText += `✅ Basis functions: ${Object.values(data.basis_functions).flat().length} types\n`;
          reportText += `✅ Exchange-correlation functionals: ${Object.values(data.exchange_correlation).flat().length} types\n`;
          reportText += `✅ Solvents: ${Object.values(data.solvents).flat().length} types\n\n`;
          
          results.dependencies.available.push('Calculation parameters');
        }
      } catch (error) {
        reportText += '❌ Parameter retrieval: Failed\n';
        reportText += `⚠️ Error: ${error instanceof Error ? error.message : 'Unknown error'}\n\n`;
        results.dependencies.missing.push('Calculation parameters');
        results.recommendations.push('Check quantum chemistry library (PySCF) initialization');
      }
    }

    // 5. StartCalculation test (only for detailed diagnostics)
    if (results.server.connected && detailed) {
      reportText += '**5. Calculation Start Function Test** (previously 500 error)\n';
      
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
          reportText += `✅ Calculation start: Normal (ID: ${calc.id})\n`;
          reportText += `   📋 Status: ${calc.status}\n`;
          
          // Try to immediately cancel the test calculation to avoid resource waste
          try {
            // Note: Cancel function not implemented in current client, so we skip this
            reportText += `   🔄 This is a test calculation, please cancel manually\n`;
          } catch (e) {
            // Ignore cancel errors
          }
          
          results.dependencies.available.push('Quantum chemistry calculation engine');
        } else {
          reportText += `❌ Calculation start: Failed\n`;
          results.dependencies.missing.push('Quantum chemistry calculation engine');
          results.recommendations.push('Issue with quantum chemistry calculation initialization');
        }
      } catch (error) {
        const errorMsg = error instanceof PySCFApiError ? 
          `${error.details.status} - ${error.details.statusText}` : 
          (error instanceof Error ? error.message : 'Unknown error');
        reportText += `❌ Calculation start: Failed (${errorMsg})\n`;
        
        if (error instanceof PySCFApiError && error.details.status === 500) {
          results.recommendations.push('Server internal error: issue with quantum chemistry library or resource management');
        }

        results.dependencies.missing.push('Quantum chemistry calculation engine');
      }
      reportText += '\n';
    }

    // 6. Summary and recommendations
    const elapsedTime = Date.now() - startTime;
    reportText += '**📊 Diagnostic Summary**\n';
    reportText += `⏱️ Diagnostic time: ${elapsedTime}ms\n`;
    reportText += `🔗 Server status: ${results.server.connected ? 'Connected' : 'Disconnected'}\n`;
    
    if (detailed) {
      reportText += `🎯 Endpoint tests: ${results.endpoints.passed}/${results.endpoints.tested} succeeded\n`;
    }
    
    reportText += `📦 Available dependencies: ${results.dependencies.available.length}\n`;
    reportText += `⚠️ Missing dependencies: ${results.dependencies.missing.length}\n\n`;

    // 7. Recommendations
    if (results.recommendations.length > 0) {
      reportText += '**🔧 Recommendations**\n';
      results.recommendations.forEach((rec, i) => {
        reportText += `${i + 1}. ${rec}\n`;
      });
      reportText += '\n';
    }

    // 8. Overall health status
    const overallHealth = results.server.connected && results.endpoints.failed === 0;
    reportText += `**🏥 Overall Health Status: ${overallHealth ? 'Healthy' : 'Needs attention'}**\n`;
    
    if (overallHealth) {
      reportText += '✅ System is operating normally.';
    } else {
      reportText += '⚠️ Some features have issues. Please check the recommendations above.';
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
          text: `❌ Error occurred during diagnostic execution: ${errorMessage}

**Problem:**
An error occurred in the diagnostic process itself.

**Solutions:**
1. Check network connection
2. Verify PySCF Native App startup status
3. Wait a moment and try again`,
        },
      ],
      isError: true,
    };
  }
}