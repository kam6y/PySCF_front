import { app, BrowserWindow, ipcMain } from 'electron'
import { join } from 'path'
import { spawn, ChildProcess } from 'child_process'
import * as fs from 'fs'

// Handle creating/removing shortcuts on Windows when installing/uninstalling.
if (require('electron-squirrel-startup')) {
  app.quit()
}

class PySCFApp {
  private mainWindow: BrowserWindow | null = null
  private pythonProcess: ChildProcess | null = null

  constructor() {
    this.setupApp()
    this.setupIPC()
  }

  private setupApp(): void {
    app.on('ready', this.createWindow.bind(this))
    
    app.on('window-all-closed', () => {
      if (process.platform !== 'darwin') {
        this.cleanup()
        app.quit()
      }
    })

    app.on('activate', () => {
      if (BrowserWindow.getAllWindows().length === 0) {
        this.createWindow()
      }
    })

    app.on('before-quit', () => {
      this.cleanup()
    })
  }

  private createWindow(): void {
    this.mainWindow = new BrowserWindow({
      height: 800,
      width: 1200,
      minHeight: 600,
      minWidth: 800,
      titleBarStyle: 'default',
      webPreferences: {
        nodeIntegration: false,
        contextIsolation: true,
        preload: join(__dirname, 'preload.js'),
      },
    })

    // Load the app
    if (process.env.NODE_ENV === 'development') {
      this.mainWindow.loadURL('http://localhost:3000')
      this.mainWindow.webContents.openDevTools()
    } else {
      this.mainWindow.loadFile(join(__dirname, '../renderer/index.html'))
    }

    this.startPythonBackend()
  }

  private startPythonBackend(): void {
    // Use absolute path to the Python executable in the project's venv
    const pythonPath = join(process.cwd(), 'venv', 'bin', 'python')
    
    const pythonScript = join(__dirname, '../python/main.py')

    // Check if Python script exists
    if (!fs.existsSync(pythonScript)) {
      console.warn('Python script not found:', pythonScript)
      return
    }

    this.pythonProcess = spawn(pythonPath, [pythonScript], {
      stdio: ['pipe', 'pipe', 'pipe'],
    })

    this.pythonProcess.stdout?.on('data', (data) => {
      console.log(`Python stdout: ${data}`)
      this.mainWindow?.webContents.send('python-log', {
        type: 'stdout',
        data: data.toString(),
      })
    })

    this.pythonProcess.stderr?.on('data', (data) => {
      console.error(`Python stderr: ${data}`)
      this.mainWindow?.webContents.send('python-log', {
        type: 'stderr',
        data: data.toString(),
      })
    })

    this.pythonProcess.on('error', (error) => {
      console.error('Failed to start Python process:', error)
      this.mainWindow?.webContents.send('python-error', error.message)
    })

    this.pythonProcess.on('close', (code) => {
      console.log(`Python process exited with code ${code}`)
      this.pythonProcess = null
    })
  }

  private setupIPC(): void {
    // Handle calculation requests
    ipcMain.handle('calculate', async (event, calculationData) => {
      return this.sendToPython('calculate', calculationData)
    })

    // Handle molecule building requests
    ipcMain.handle('build-molecule', async (event, moleculeData) => {
      return this.sendToPython('build-molecule', moleculeData)
    })

    // Handle getting calculation status
    ipcMain.handle('get-calculation-status', async () => {
      return this.sendToPython('get-status', {})
    })

    // Handle stopping calculations
    ipcMain.handle('stop-calculation', async () => {
      return this.sendToPython('stop', {})
    })
  }

  private sendToPython(action: string, data: any): Promise<any> {
    return new Promise((resolve, reject) => {
      if (!this.pythonProcess) {
        reject(new Error('Python process not available'))
        return
      }

      const message = JSON.stringify({ action, data, id: Date.now() })
      
      // Set up one-time listener for response
      const responseHandler = (response: Buffer) => {
        try {
          const result = JSON.parse(response.toString())
          resolve(result)
        } catch (error) {
          reject(error)
        }
      }

      this.pythonProcess.stdout?.once('data', responseHandler)
      this.pythonProcess.stdin?.write(message + '\n')

      // Timeout after 30 seconds
      setTimeout(() => {
        this.pythonProcess?.stdout?.removeListener('data', responseHandler)
        reject(new Error('Python process timeout'))
      }, 30000)
    })
  }

  private cleanup(): void {
    if (this.pythonProcess) {
      this.pythonProcess.kill()
      this.pythonProcess = null
    }
  }
}

// Initialize the application
new PySCFApp()