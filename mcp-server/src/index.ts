#!/usr/bin/env node

import { Server } from '@modelcontextprotocol/sdk/server/index.js';
import { StdioServerTransport } from '@modelcontextprotocol/sdk/server/stdio.js';
import {
  CallToolRequestSchema,
  ListToolsRequestSchema,
  ErrorCode,
  McpError,
} from '@modelcontextprotocol/sdk/types.js';
import { PySCFApiClient } from './client.js';
import { tools, toolHandlers } from './tools/index.js';
import { MCPInspector } from './tools/inspector.js';

class PySCFMCPServer {
  private server: Server;
  private client: PySCFApiClient | null = null;

  constructor() {
    this.server = new Server(
      {
        name: 'pyscf-mcp-server',
        version: '1.0.0',
      },
      {
        capabilities: {
          tools: {},
        },
      }
    );

    this.setupToolHandlers();
    this.setupErrorHandler();
  }

  private setupErrorHandler(): void {
    this.server.onerror = (error) => {
      console.error('[MCP Error]', error);
    };

    process.on('SIGINT', async () => {
      await this.server.close();
      process.exit(0);
    });
  }

  private setupToolHandlers(): void {
    this.server.setRequestHandler(ListToolsRequestSchema, async () => {
      return {
        tools: tools,
      };
    });

    this.server.setRequestHandler(CallToolRequestSchema, async (request) => {
      const { name, arguments: args } = request.params;
      const startTime = Date.now();

      try {
        // Ensure we have a connected client
        if (!this.client) {
          await this.initializeClient();
        }

        // Check if tool exists
        const handler = toolHandlers[name];
        if (!handler) {
          throw new McpError(
            ErrorCode.MethodNotFound,
            `Unknown tool: ${name}`
          );
        }

        // Execute the tool
        const result = await handler(args || {}, this.client!);
        
        // Log successful execution
        MCPInspector.addLog({
          timestamp: new Date(),
          tool: name,
          args: args || {},
          success: true,
          responseTime: Date.now() - startTime,
        });
        
        return result;

      } catch (error) {
        console.error(`Error executing tool ${name}:`, error);
        
        // Log failed execution
        const errorMessage = error instanceof Error ? error.message : 'Unknown error';
        MCPInspector.addLog({
          timestamp: new Date(),
          tool: name,
          args: args || {},
          success: false,
          error: errorMessage,
          responseTime: Date.now() - startTime,
        });

        if (error instanceof McpError) {
          throw error;
        }

        // Handle connection errors specifically
        if (error instanceof Error && error.message.includes('ECONNREFUSED')) {
          throw new McpError(
            ErrorCode.InternalError,
            `PySCF Native Appサーバーに接続できません。アプリが起動していることを確認してください。\n\nエラー: ${error.message}`
          );
        }

        throw new McpError(
          ErrorCode.InternalError,
          `Tool execution failed: ${error instanceof Error ? error.message : 'Unknown error'}`
        );
      }
    });
  }

  private async initializeClient(): Promise<void> {
    try {
      console.error('PySCF Native Appサーバーに接続中...');
      
      // Try auto-detection first
      this.client = await PySCFApiClient.createWithAutoDetect('127.0.0.1', {
        start: 5000,
        end: 5100,
      });

      const connection = await this.client.testConnection();
      if (connection.connected) {
        console.error(`✅ PySCF Native Appに接続成功: ${this.client.getBaseUrl()}`);
        console.error(`   サーバーバージョン: ${connection.version || 'N/A'}`);
      } else {
        throw new Error(connection.error || 'Connection test failed');
      }

    } catch (error) {
      console.error('❌ PySCF Native Appサーバーへの接続に失敗しました');
      console.error('');
      console.error('解決方法:');
      console.error('1. PySCF Native Appが起動していることを確認してください');
      console.error('2. ターミナルで以下のコマンドを実行してアプリを起動:');
      console.error('   cd /path/to/PySCF_native_app');
      console.error('   npm run dev');
      console.error('3. アプリが http://127.0.0.1:5000 で起動していることを確認');
      console.error('');
      
      throw new Error(
        `PySCF Native Appサーバーが見つかりません。\n` +
        `アプリを起動してから再試行してください。\n\n` +
        `詳細: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  async run(): Promise<void> {
    const transport = new StdioServerTransport();
    
    try {
      console.error('PySCF MCP Server starting...');
      console.error('Available tools:', tools.length);
      console.error('Tool names:', tools.map(t => t.name).join(', '));
      
      await this.server.connect(transport);
      console.error('✅ MCP Server ready');
      
    } catch (error) {
      console.error('Failed to start MCP server:', error);
      process.exit(1);
    }
  }
}

// Start the server
const server = new PySCFMCPServer();
server.run().catch((error) => {
  console.error('Fatal error:', error);
  process.exit(1);
});