import path from "node:path";
import { BrowserWindow, app } from "electron";
import { spawn, ChildProcess } from "child_process";
import http from "node:http";

let mainWindow: BrowserWindow;
let pythonProcess: ChildProcess | null = null;

/**
 * Pythonサーバーのヘルスチェックを行い、起動完了を待つ関数
 * @param url チェック対象のURL
 * @param retries リトライ回数
 * @param delay リトライ間隔（ミリ秒）
 * @returns サーバーが正常に応答すれば解決されるPromise
 */
const checkServerHealth = (url: string, retries = 15, delay = 1000): Promise<void> => {
  return new Promise((resolve, reject) => {
    let attempts = 0;
    const interval = setInterval(() => {
      http.get(url, (res) => {
        if (res.statusCode === 200) {
          clearInterval(interval);
          resolve();
        }
      }).on('error', () => {
        attempts++;
        if (attempts >= retries) {
          clearInterval(interval);
          reject(new Error('Python server health check failed after multiple retries.'));
        }
      });
    }, delay);
  });
};

/**
 * Python/Flaskサーバーを起動する関数
 * app.isPackaged の値に応じて、開発モードとパッケージ後で起動方法を切り替える
 */
const startPythonServer = async (): Promise<void> => {
  console.log('Starting Python Flask server...');
  
  if (app.isPackaged) {
    // パッケージ化後の場合：同梱した実行ファイルを起動
    const executableName = process.platform === 'win32' ? 'pyscf_api.exe' : 'pyscf_api';
    // extraResources は process.resourcesPath にコピーされる
    const pythonExecutablePath = path.join(process.resourcesPath, 'python_dist', executableName);
    
    console.log(`Executing packaged Python server at: ${pythonExecutablePath}`);
    pythonProcess = spawn(pythonExecutablePath, [], {
      stdio: ['pipe', 'pipe', 'pipe']
    });

  } else {
    // 開発時の場合：uv を使って直接Pythonスクリプトを起動
    const pythonPath = path.join(__dirname, '..', 'src', 'python');
    const appPath = path.join(pythonPath, 'app.py');
    
    console.log('Executing Python server in development mode...');
    pythonProcess = spawn('uv', ['run', 'python', appPath], {
      cwd: pythonPath,
      stdio: ['pipe', 'pipe', 'pipe']
    });
  }

  // stdout/stderrのログ出力
  pythonProcess.stdout?.on('data', (data) => {
    console.log(`Python server: ${data.toString().trim()}`);
  });

  pythonProcess.stderr?.on('data', (data) => {
    console.error(`Python server error: ${data.toString().trim()}`);
  });

  // プロセス自体のエラーハンドリング
  pythonProcess.on('error', (error) => {
    console.error(`Failed to start Python server process: ${error.message}`);
  });

  pythonProcess.on('close', (code) => {
    console.log(`Python server exited with code ${code}`);
    pythonProcess = null;
  });

  // サーバーのヘルスチェックを行い、起動を待つ
  await checkServerHealth('http://127.0.0.1:5000/health');
};

/**
 * Python/Flaskサーバーを停止する関数
 */
const stopPythonServer = (): void => {
  if (pythonProcess && !pythonProcess.killed) {
    console.log('Stopping Python Flask server...');
    // まずは穏便に終了を試みる (SIGTERM)
    pythonProcess.kill('SIGTERM'); 
    
    // 5秒経っても終了しない場合は強制終了 (SIGKILL)
    setTimeout(() => {
      if (pythonProcess && !pythonProcess.killed) {
        console.log('Force killing Python server...');
        pythonProcess.kill('SIGKILL');
      }
    }, 5000);
  }
};


// Electronアプリの準備が完了したら実行
app.whenReady().then(async () => {
  // 最初にPythonサーバーを起動
  try {
    await startPythonServer();
    console.log('Python server started successfully and is healthy.');
  } catch (error) {
    console.error('Failed to start Python server:', error);
    // サーバーが起動しない場合はアプリを終了
    app.quit();
    return;
  }
  
  // BrowserWindow インスタンスを作成
  mainWindow = new BrowserWindow({
    width: 1200,
    height: 800,
    minWidth: 800,
    minHeight: 600,
    titleBarStyle: 'hidden',
    titleBarOverlay: {
      color: 'rgba(0, 0, 0, 0)',
      symbolColor: '#000000',
      height: 40
    },
    webPreferences: {
      preload: path.join(__dirname, "preload.js"),
      nodeIntegration: false,
      contextIsolation: true,
    },
  });

  // 開発モードの場合はデベロッパーツールを開く
  if (!app.isPackaged) {
    mainWindow.webContents.openDevTools({ mode: 'detach' });
  }

  // レンダラープロセスをロード
  mainWindow.loadFile("dist/index.html");
});


// すべてのウィンドウが閉じられたらアプリを終了する
app.once("window-all-closed", () => {
  stopPythonServer();
  app.quit();
});


// アプリ終了前の処理
app.on('before-quit', () => {
  stopPythonServer();
});