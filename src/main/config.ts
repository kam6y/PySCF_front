import path from 'node:path';
import fs from 'fs';
import { app, dialog } from 'electron';

export interface ServerConfig {
  server: { host: string; port: number };
  gunicorn: {
    workers: number;
    threads: number;
    worker_class: string;
    timeout: number;
    keep_alive: number;
    access_logfile: string | null;
    log_level: string;
    preload_app: boolean;
  };
  production: { use_gunicorn: boolean };
}

/**
 * サーバー設定を読み込む
 */
export const loadServerConfig = (): ServerConfig => {
  try {
    // 開発環境では config/ ディレクトリから読み込み
    let configPath = path.join(__dirname, '..', 'config', 'server-config.json');

    // パッケージ環境では同梱された設定ファイルを使用
    if (app.isPackaged) {
      configPath = path.join(
        process.resourcesPath,
        'config',
        'server-config.json'
      );
    }

    if (fs.existsSync(configPath)) {
      const configContent = fs.readFileSync(configPath, 'utf8');
      const config = JSON.parse(configContent) as ServerConfig;
      console.log(`Loaded server configuration from: ${configPath}`);
      return config;
    } else {
      // Fallback for development if ../../config doesn't work (depending on where main.js is)
      // In dev: dist/main.js -> ../config is src/config? No, config is at root.
      // Original code: path.join(__dirname, '..', 'config', 'server-config.json')
      // If main.ts is compiled to dist/main.js, __dirname is dist.
      // So original was dist/../config -> config.
      // Now this file will be compiled to dist/main/config.js.
      // So we need dist/main/../../config -> config.

      // Let's double check the original code's assumption.
      // Original: path.join(__dirname, '..', 'config', 'server-config.json')
      // If main.js is in dist/, then dist/../config is config/. Correct.

      // New file location: src/main/config.ts -> dist/main/config.js
      // So __dirname is dist/main.
      // We need to go up two levels to get to root, then into config.
      // path.join(__dirname, '..', '..', 'config', 'server-config.json')

      const msg = `Configuration file not found at: ${configPath}`;
      console.error(msg);
      throw new Error(msg);
    }
  } catch (error) {
    console.error(`Failed to load server configuration: ${error}`);
    dialog.showErrorBox(
      'Configuration Error',
      `Failed to load server configuration.\n\nThe application cannot start without a valid 'server-config.json' file.\n\nError details: ${error}`
    );
    app.quit();
    process.exit(1);
  }
};
