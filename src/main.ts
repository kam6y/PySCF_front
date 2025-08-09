import path from "node:path";
import { BrowserWindow, app } from "electron";
import { spawn, ChildProcess } from "child_process";

let mainWindow: BrowserWindow;
let pythonProcess: ChildProcess | null = null;

const startPythonServer = async (): Promise<void> => {
  return new Promise((resolve, reject) => {
    console.log('Starting Python Flask server...');
    
    const pythonPath = path.join(__dirname, '..', 'src', 'python');
    const appPath = path.join(pythonPath, 'app.py');
    
    // Use uv to run the Python server
    pythonProcess = spawn('uv', ['run', 'python', appPath], {
      cwd: pythonPath,
      stdio: ['pipe', 'pipe', 'pipe']
    });

    pythonProcess.stdout?.on('data', (data) => {
      console.log(`Python server: ${data.toString()}`);
      // Look for server startup message
      if (data.toString().includes('Running on')) {
        resolve();
      }
    });

    pythonProcess.stderr?.on('data', (data) => {
      console.error(`Python server error: ${data.toString()}`);
    });

    pythonProcess.on('error', (error) => {
      console.error(`Failed to start Python server: ${error.message}`);
      reject(error);
    });

    pythonProcess.on('close', (code) => {
      console.log(`Python server exited with code ${code}`);
      pythonProcess = null;
    });

    // Timeout after 10 seconds
    setTimeout(() => {
      if (pythonProcess && !pythonProcess.killed) {
        resolve(); // Assume it started even if we didn't see the message
      }
    }, 10000);
  });
};

const stopPythonServer = (): void => {
  if (pythonProcess && !pythonProcess.killed) {
    console.log('Stopping Python Flask server...');
    pythonProcess.kill('SIGTERM');
    
    // Force kill after 5 seconds if still running
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
    console.log('Python server started successfully');
  } catch (error) {
    console.error('Failed to start Python server:', error);
    // Continue anyway - maybe the server is already running
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