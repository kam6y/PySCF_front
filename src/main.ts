import path from "node:path";
import { BrowserWindow, app } from "electron";
import { spawn, ChildProcess } from "child_process";
import http from "node:http";

let mainWindow: BrowserWindow;
let pythonProcess: ChildProcess | null = null;

//
// 変更点 1: パス解決の堅牢化
// app.isPackaged を使用して、開発時とパッケージ化後でパスを切り替える
//
const getPythonPath = () => {
  if (app.isPackaged) {
    // パッケージ化後は、リソースディレクトリに配置されたPythonフォルダを参照する
    // electron-builderのextraResources設定で 'src/python' をコピーする必要がある
    return path.join(process.resourcesPath, 'python');
  } else {
    // 開発時は従来のパスを使用
    return path.join(__dirname, '..', 'src', 'python');
  }
};

//
// 変更点 2: ヘルスチェックによる起動確認
//
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


const startPythonServer = async (): Promise<void> => {
  console.log('Starting Python Flask server...');
  
  const pythonPath = getPythonPath();
  const appPath = path.join(pythonPath, 'app.py');
  
  // Use uv to run the Python server
  pythonProcess = spawn('uv', ['run', 'python', appPath], {
    cwd: pythonPath,
    stdio: ['pipe', 'pipe', 'pipe'] // stdioをパイプしてログを確認可能にする
  });

  pythonProcess.stdout?.on('data', (data) => {
    console.log(`Python server: ${data.toString().trim()}`);
  });

  pythonProcess.stderr?.on('data', (data) => {
    console.error(`Python server error: ${data.toString().trim()}`);
  });

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

const stopPythonServer = (): void => {
  if (pythonProcess && !pythonProcess.killed) {
    console.log('Stopping Python Flask server...');
    pythonProcess.kill('SIGTERM'); // まずは穏便に終了を試みる
    
    // 5秒経っても終了しない場合は強制終了
    setTimeout(() => {
      if (pythonProcess && !pythonProcess.killed) {
        console.log('Force killing Python server...');
        pythonProcess.kill('SIGKILL');
      }
    }, 5000);
  }
};

app.whenReady().then(async () => {
  // Start Python server first
  try {
    await startPythonServer();
    console.log('Python server started successfully and is healthy.');
  } catch (error) {
    console.error('Failed to start Python server:', error);
    // サーバーが起動しない場合、アプリを終了するか、エラーダイアログを表示するなどの処理が考えられる
    // ここでは、とりあえず続行せずにアプリを終了させる
    app.quit();
    return;
  }
  // アプリの起動イベント発火で BrowserWindow インスタンスを作成
  mainWindow = new BrowserWindow({
    width: 1200,
    height: 800,
    minWidth: 800,
    minHeight: 600,
    titleBarStyle: 'hidden', // Hide native title bar
    titleBarOverlay: {
      color: 'rgba(0, 0, 0, 0)', // Transparent overlay
      symbolColor: '#000000',
      height: 40
    },
    webPreferences: {
      // webpack が出力したプリロードスクリプトを読み込み
      preload: path.join(__dirname, "preload.js"),
      nodeIntegration: false,
      contextIsolation: true,
    },
  });

  mainWindow.webContents.openDevTools({ mode: 'detach' });
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