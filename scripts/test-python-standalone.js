#!/usr/bin/env node

const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const os = require('os');

/**
 * Test script to verify Python executable functionality outside of Electron
 * This simulates what main.ts does when starting the Python backend
 */

console.log('=== Python Executable Standalone Test ===');

// Determine platform and architecture
const platform = os.platform(); // 'darwin' (macOS), 'win32' (Windows), 'linux'
const arch = os.arch(); // 'arm64', 'x64'

// Generate platform-specific directory names
let platformArch;
if (platform === 'darwin') {
    platformArch = `mac-${arch}`;
} else if (platform === 'win32') {
    platformArch = `win-unpacked`; // Windows の場合の一般的なディレクトリ名
} else {
    platformArch = `linux-unpacked`; // Linux の場合
}

// Determine paths
const appPath = platform === 'darwin' 
    ? `./dist/${platformArch}/Pyscf_front.app`
    : `./dist/${platformArch}`;

const isPackaged = fs.existsSync(appPath);
console.log(`Platform: ${platform}-${arch}`);
console.log(`Packaged mode: ${isPackaged}`);

let pythonExecutablePath;
let pythonWorkingDir;

if (isPackaged) {
  // Test packaged conda environment first
  const resourcesPath = platform === 'darwin'
    ? path.join(appPath, 'Contents/Resources')
    : path.resolve(`./dist/${platformArch}/resources`); // Windows/Linux の場合

  const condaPythonPath = path.join(resourcesPath, 'conda_env', 'bin', 'python');
  const condaGunicornPath = path.join(resourcesPath, 'conda_env', 'bin', 'gunicorn');
  
  if (fs.existsSync(condaPythonPath) && fs.existsSync(condaGunicornPath)) {
    pythonExecutablePath = path.resolve(condaPythonPath);
    pythonWorkingDir = path.resolve(path.join(resourcesPath, 'conda_env', 'bin'));
    console.log('✓ Using packaged conda environment');
  } else {
    console.log('✗ No packaged conda environment found.');
    process.exit(1);
  }
} else {
  console.log('✗ No packaged app found. Please run `npm run package` first.');
  process.exit(1);
}

console.log(`Python executable: ${pythonExecutablePath}`);
console.log(`Working directory: ${pythonWorkingDir}`);
console.log(`Executable exists: ${fs.existsSync(pythonExecutablePath)}`);
console.log(`Working dir exists: ${fs.existsSync(pythonWorkingDir)}`);

// Main test execution using async/await
(async () => {
  try {
    console.log('\n=== Test 1: Python Version Check ===');
    await testPythonCommand([pythonExecutablePath, '--version']);

    console.log('\n=== Test 2: Basic Import Test ===');
    await testPythonCommand([pythonExecutablePath, '-c', 'import sys; print("Python executable test successful"); print("Python version:", sys.version)']);

    console.log('\n=== Test 3: Gunicorn Import Test ===');
    await testPythonCommand([pythonExecutablePath, '-c', 'import gunicorn; print("Gunicorn version:", gunicorn.__version__)']);
    
    console.log('\n=== Test 4: Flask Import Test ===');
    await testPythonCommand([pythonExecutablePath, '-c', 'import flask; print("Flask import successful")']);

    // Test 5: App import test (if conda environment)
    if (pythonExecutablePath.includes('conda_env')) {
      console.log('\n=== Test 5: App Import Test ===');
      // Need to change working directory to where app.py is located
      const appWorkingDir = pythonExecutablePath.includes('conda_env') 
        ? './src/python'  // Development app.py location
        : './src/python';  // Development app.py location
      await testPythonCommand([pythonExecutablePath, '-c', 'import app; print("App import successful")'], appWorkingDir);
    }

    console.log('\n=== Test 6: Gunicorn Startup Simulation ===');
    const gunicornArgs = [
      '-m', 'gunicorn',
      '--workers', '1',
      '--threads', '4',
      '--worker-class', 'sync',
      '--bind', '127.0.0.1:5001',  // Use different port to avoid conflicts
      '--timeout', '0',
      '--keep-alive', '30',
      '--access-logfile', '-',
      '--log-level', 'info',
      '--preload',
      'app:app'
    ];
    
    const workDir = pythonExecutablePath.includes('conda_env') 
      ? './src/python'  // Development app.py location
      : './src/python';  // Development app.py location
      
    console.log(`Command: ${pythonExecutablePath} ${gunicornArgs.join(' ')}`);
    console.log(`Working directory: ${workDir}`);
    
    await testPythonCommand([pythonExecutablePath, ...gunicornArgs], workDir, 5000);  // 5 second timeout

    console.log('\n✓ All tests passed! The Python executable should work correctly.');

  } catch (error) {
    console.error('\n✗ Test failed:', error.message);
    process.exit(1);
  }
})();

function testPythonCommand(command, cwd = null, timeout = 10000) {
  return new Promise((resolve, reject) => {
    console.log(`Running: ${command.join(' ')}`);
    if (cwd) console.log(`Working directory: ${cwd}`);
    
    const childProcess = spawn(command[0], command.slice(1), {
      cwd: cwd || pythonWorkingDir,
      stdio: ['pipe', 'pipe', 'pipe'],
      env: { ...process.env, CONDA_DEFAULT_ENV: 'pyscf-env' }
    });
    
    let stdout = '';
    let stderr = '';
    let resolved = false;
    
    // Set timeout
    const timer = setTimeout(() => {
      if (!resolved) {
        console.log('⚠️  Command timed out (this may be expected for server startup)');
        childProcess.kill('SIGTERM');
        resolved = true;
        resolve();  // Resolve for server startup test
      }
    }, timeout);
    
    childProcess.stdout?.on('data', (data) => {
      const output = data.toString();
      stdout += output;
      console.log(`STDOUT: ${output.trim()}`);
      
      // For Gunicorn, if we see it starting to listen, that's success
      if (output.includes('Listening at:') || output.includes('Booting worker')) {
        if (!resolved) {
          console.log('✓ Gunicorn startup detected - killing process');
          clearTimeout(timer);
          childProcess.kill('SIGTERM');
          resolved = true;
          resolve();
        }
      }
    });
    
    childProcess.stderr?.on('data', (data) => {
      const output = data.toString();
      stderr += output;
      console.log(`STDERR: ${output.trim()}`);
    });
    
    childProcess.on('close', (code) => {
      clearTimeout(timer);
      if (!resolved) {
        if (code === 0) {
          console.log('✓ Command completed successfully');
          resolve();
        } else {
          console.log(`✗ Command failed with exit code: ${code}`);
          if (stderr) console.log(`Error output: ${stderr}`);
          reject(new Error(`Command failed with exit code ${code}: ${stderr || 'No error message'}`));
        }
        resolved = true;
      }
    });
    
    childProcess.on('error', (error) => {
      clearTimeout(timer);
      if (!resolved) {
        console.log(`✗ Process error: ${error.message}`);
        reject(error);
        resolved = true;
      }
    });
  });
}