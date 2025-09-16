#!/usr/bin/env node

/**
 * MCP設定ファイル生成スクリプト
 * 
 * このスクリプトは、Claude Desktop用のMCP設定ファイルを
 * テンプレートから自動生成します。プロジェクトパスを自動検出し、
 * ハードコードされたパスの問題を解決します。
 */

const fs = require('fs');
const path = require('path');
const process = require('process');

// カラー出力用のヘルパー
const colors = {
  green: '\x1b[32m',
  yellow: '\x1b[33m',
  red: '\x1b[31m',
  blue: '\x1b[34m',
  reset: '\x1b[0m',
  bold: '\x1b[1m'
};

function log(message, color = 'reset') {
  console.log(`${colors[color]}${message}${colors.reset}`);
}

function logSuccess(message) {
  log(`✓ ${message}`, 'green');
}

function logInfo(message) {
  log(`ℹ ${message}`, 'blue');
}

function logWarning(message) {
  log(`⚠ ${message}`, 'yellow');
}

function logError(message) {
  log(`✗ ${message}`, 'red');
}

function findProjectRoot() {
  let currentDir = __dirname;
  
  // scriptsディレクトリから親ディレクトリに移動
  const projectRoot = path.dirname(currentDir);
  
  // package.jsonの存在確認
  const packageJsonPath = path.join(projectRoot, 'package.json');
  if (!fs.existsSync(packageJsonPath)) {
    throw new Error('プロジェクトルートが見つかりません (package.jsonが見つかりません)');
  }
  
  // mcp-serverディレクトリの存在確認  
  const mcpServerPath = path.join(projectRoot, 'mcp-server');
  if (!fs.existsSync(mcpServerPath)) {
    throw new Error('mcp-serverディレクトリが見つかりません');
  }
  
  return projectRoot;
}

function generateConfig() {
  try {
    logInfo('MCP設定ファイル生成を開始します...');
    
    // プロジェクトルートを検出
    const projectRoot = findProjectRoot();
    logInfo(`プロジェクトルート: ${projectRoot}`);
    
    // テンプレートファイルパス
    const templatePath = path.join(projectRoot, 'mcp-server', 'claude_desktop_config.template.json');
    const outputPath = path.join(projectRoot, 'mcp-server', 'claude_desktop_config.json');
    
    // テンプレートファイルの存在確認
    if (!fs.existsSync(templatePath)) {
      throw new Error(`テンプレートファイルが見つかりません: ${templatePath}`);
    }
    
    // テンプレートファイルを読み込み
    const templateContent = fs.readFileSync(templatePath, 'utf8');
    
    // プレースホルダーを実際のパスに置換
    const configContent = templateContent.replace(/\{\{PROJECT_PATH\}\}/g, projectRoot);
    
    // 設定ファイルを出力
    fs.writeFileSync(outputPath, configContent, 'utf8');
    
    logSuccess('設定ファイルが正常に生成されました');
    logInfo(`出力先: ${outputPath}`);
    
    // 生成された設定ファイルの内容を表示
    logInfo('\\n生成された設定内容:');
    log(configContent, 'blue');
    
    // Claude Desktopの設定手順を表示
    logInfo('\\n次のステップ:');
    log('1. Claude Desktop設定ファイルに以下の内容を追加してください:', 'yellow');
    log('   場所: ~/Library/Application Support/Claude/claude_desktop_config.json', 'yellow');
    log('2. 上記の生成された設定内容をコピーして設定ファイルに貼り付けてください', 'yellow');
    log('3. Claude Desktopを再起動してください', 'yellow');
    
    logSuccess('\\n設定ファイル生成が完了しました!');
    
  } catch (error) {
    logError(`エラーが発生しました: ${error.message}`);
    process.exit(1);
  }
}

// ヘルプ表示
function showHelp() {
  log('\\nMCP設定ファイル生成スクリプト', 'bold');
  log('=====================================\\n');
  
  log('使用方法:', 'blue');
  log('  npm run setup-mcp-config    # 設定ファイルを生成');
  log('  node scripts/setup-mcp-config.js --help    # ヘルプを表示\\n');
  
  log('このスクリプトは以下を実行します:', 'blue');
  log('• プロジェクトパスを自動検出');
  log('• テンプレートファイルからClaude Desktop設定を生成');
  log('• ハードコードされたパスの問題を解決\\n');
}

// コマンドライン引数の処理
if (process.argv.includes('--help') || process.argv.includes('-h')) {
  showHelp();
  process.exit(0);
}

// メイン処理実行
generateConfig();