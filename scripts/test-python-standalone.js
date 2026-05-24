#!/usr/bin/env node

const { spawn, spawnSync } = require('child_process');
const path = require('path');
const fs = require('fs');
const os = require('os');
const http = require('http');
const net = require('net');
const packageJson = require('../package.json');

/**
 * Test script to verify Python executable functionality outside of Electron
 * This simulates what main.ts does when starting the Python backend
 */

console.log('=== Python Executable Standalone Test ===');

const AUTH_TOKEN = 'standalone-test-token';
const SERVER_ENV = 'production';
const CONDA_ENV_MARKER_FILE = '.pyscf-standalone-conda-env';

// Determine platform and architecture
const platform = os.platform(); // 'darwin' (macOS), 'win32' (Windows), 'linux'
const arch = os.arch(); // 'arm64', 'x64'

function getPlatformDirectories() {
  if (platform === 'darwin') {
    return [`mac-${arch}`, 'mac'];
  }
  if (platform === 'win32') {
    return ['win-unpacked'];
  }
  return ['linux-unpacked'];
}

function getAppPath(outputRoot, platformDirectory) {
  if (platform === 'darwin') {
    return path.resolve(outputRoot, platformDirectory, 'Pyscf_front.app');
  }
  return path.resolve(outputRoot, platformDirectory);
}

function getResourcesPath(appPath) {
  if (platform === 'darwin') {
    return path.resolve(appPath, 'Contents/Resources');
  }
  return path.resolve(appPath, 'resources');
}

function findPackagedApp() {
  const outputRoots = ['release', 'dist'];
  const platformDirectories = getPlatformDirectories();

  for (const outputRoot of outputRoots) {
    for (const platformDirectory of platformDirectories) {
      const appPath = getAppPath(outputRoot, platformDirectory);
      if (fs.existsSync(appPath)) {
        return {
          appPath,
          outputRoot,
          platformDirectory,
          resourcesPath: getResourcesPath(appPath),
        };
      }
    }
  }

  return null;
}

const packagedApp = findPackagedApp();
const isPackaged = packagedApp !== null;
console.log(`Platform: ${platform}-${arch}`);
console.log(`Packaged mode: ${isPackaged}`);

let pythonExecutablePath;
let pythonWorkingDir;
let resourcesPath;
let packagedPythonSourceDir;
let runtimeCondaDir;

if (isPackaged) {
  const { outputRoot, platformDirectory } = packagedApp;
  console.log(`Packaged output root: ${outputRoot}`);
  console.log(`Packaged directory: ${platformDirectory}`);

  // Test packaged conda environment first
  resourcesPath = packagedApp.resourcesPath;

  const bundledCondaDir = path.join(resourcesPath, 'conda_env');
  runtimeCondaDir = prepareRelocatedCondaEnvironment(bundledCondaDir);
  const condaPythonPath = path.join(runtimeCondaDir, 'bin', 'python');
  const condaGunicornPath = path.join(runtimeCondaDir, 'bin', 'gunicorn');
  packagedPythonSourceDir = path.join(resourcesPath, 'src', 'python');

  if (fs.existsSync(condaPythonPath) && fs.existsSync(condaGunicornPath)) {
    pythonExecutablePath = path.resolve(condaPythonPath);
    pythonWorkingDir = path.resolve(path.join(runtimeCondaDir, 'bin'));
    console.log('✓ Using relocated packaged conda environment');
  } else {
    console.log('✗ No relocated packaged conda environment found.');
    process.exit(1);
  }
} else {
  console.log('✗ No packaged app found. Please run `npm run package` first.');
  process.exit(1);
}

console.log(`Python executable: ${pythonExecutablePath}`);
console.log(`Working directory: ${pythonWorkingDir}`);
console.log(`Packaged resources: ${resourcesPath}`);
console.log(`Packaged Python source: ${packagedPythonSourceDir}`);
console.log(`Runtime conda environment: ${runtimeCondaDir}`);
console.log(`Executable exists: ${fs.existsSync(pythonExecutablePath)}`);
console.log(`Working dir exists: ${fs.existsSync(pythonWorkingDir)}`);
console.log(`Python source exists: ${fs.existsSync(packagedPythonSourceDir)}`);

if (!fs.existsSync(packagedPythonSourceDir)) {
  console.log('✗ Packaged Python source directory not found.');
  process.exit(1);
}

function condaEnvExists(condaDir) {
  return (
    fs.existsSync(path.join(condaDir, 'bin', 'python')) &&
    fs.existsSync(path.join(condaDir, 'bin', 'gunicorn'))
  );
}

function getCondaEnvMarker(sourceCondaDir) {
  const historyPath = path.join(sourceCondaDir, 'conda-meta', 'history');
  const markerSource = fs.existsSync(historyPath)
    ? fs.statSync(historyPath)
    : fs.statSync(path.join(sourceCondaDir, 'bin', 'python'));
  return [
    packageJson.version,
    platform,
    arch,
    markerSource.size,
    Math.trunc(markerSource.mtimeMs),
  ].join(':');
}

function runCondaUnpack(condaDir) {
  const condaUnpackPath = path.join(condaDir, 'bin', 'conda-unpack');
  if (!fs.existsSync(condaUnpackPath)) {
    console.log('conda-unpack not found; using copied conda environment as-is');
    return;
  }

  console.log(`Running conda-unpack in relocated environment: ${condaDir}`);
  const condaBinDir = path.join(condaDir, 'bin');
  const result = spawnSync(condaUnpackPath, [], {
    cwd: condaDir,
    encoding: 'utf8',
    env: {
      ...process.env,
      PATH: `${condaBinDir}:${process.env.PATH || ''}`,
      CONDA_DEFAULT_ENV: 'pyscf-env',
    },
  });

  if (result.stdout)
    console.log(`conda-unpack STDOUT: ${result.stdout.trim()}`);
  if (result.stderr)
    console.log(`conda-unpack STDERR: ${result.stderr.trim()}`);
  if (result.status !== 0) {
    throw new Error(`conda-unpack failed with exit code ${result.status}`);
  }
}

function prepareRelocatedCondaEnvironment(bundledCondaDir) {
  if (!condaEnvExists(bundledCondaDir)) {
    throw new Error(
      `Bundled conda environment is missing or incomplete: ${bundledCondaDir}`
    );
  }

  const runtimeBase = path.join(
    os.tmpdir(),
    'pyscf-front-standalone',
    packagedApp.platformDirectory
  );
  const runtimeCondaPath = path.join(runtimeBase, 'conda_env');
  const expectedMarker = getCondaEnvMarker(bundledCondaDir);
  const markerPath = path.join(runtimeCondaPath, CONDA_ENV_MARKER_FILE);
  const currentMarker = fs.existsSync(markerPath)
    ? fs.readFileSync(markerPath, 'utf8').trim()
    : null;

  if (currentMarker === expectedMarker && condaEnvExists(runtimeCondaPath)) {
    console.log(`✓ Reusing relocated conda environment: ${runtimeCondaPath}`);
    return runtimeCondaPath;
  }

  console.log(
    `Copying conda environment to runtime location: ${runtimeCondaPath}`
  );
  fs.rmSync(runtimeCondaPath, { recursive: true, force: true });
  fs.mkdirSync(runtimeBase, { recursive: true });
  fs.cpSync(bundledCondaDir, runtimeCondaPath, {
    recursive: true,
    force: true,
  });
  runCondaUnpack(runtimeCondaPath);
  fs.writeFileSync(markerPath, `${expectedMarker}\n`);
  return runtimeCondaPath;
}

function createPythonEnv(serverPort = 0) {
  const condaBinDir = path.dirname(pythonExecutablePath);
  return {
    ...process.env,
    PATH: `${condaBinDir}:${process.env.PATH || ''}`,
    CONDA_DEFAULT_ENV: 'pyscf-env',
    PYSCF_RESOURCES_PATH: resourcesPath,
    PYSCF_AUTH_TOKEN: AUTH_TOKEN,
    PYSCF_ENV: SERVER_ENV,
    PYSCF_SERVER_PORT: String(serverPort),
    PYSCF_PARENT_PID: String(process.pid),
    ...(platform === 'darwin' && arch === 'arm64'
      ? { OPENBLAS_CORETYPE: 'ARMV8', LC_ALL: process.env.LC_ALL || 'C' }
      : {}),
  };
}

// Main test execution using async/await
(async () => {
  try {
    console.log('\n=== Test 1: Python Version Check ===');
    await testPythonCommand([pythonExecutablePath, '--version']);

    console.log('\n=== Test 2: Basic Import Test ===');
    await testPythonCommand([
      pythonExecutablePath,
      '-c',
      'import sys; print("Python executable test successful"); print("Python version:", sys.version)',
    ]);

    console.log('\n=== Test 3: Gunicorn Import Test ===');
    await testPythonCommand([
      pythonExecutablePath,
      '-c',
      'import gunicorn; print("Gunicorn version:", gunicorn.__version__)',
    ]);

    console.log('\n=== Test 4: FastAPI ASGI Import Test ===');
    await testPythonCommand([
      pythonExecutablePath,
      '-c',
      'import fastapi, uvicorn, socketio; from uvicorn.workers import UvicornWorker; print("FastAPI ASGI import successful")',
    ]);

    // Test 5: App import test (if conda environment)
    if (pythonExecutablePath.includes('conda_env')) {
      console.log('\n=== Test 5: App Import Test ===');
      // Need to change working directory to where packaged app.py is located
      const appWorkingDir = packagedPythonSourceDir;
      await testPythonCommand(
        [
          pythonExecutablePath,
          '-c',
          'import app; print("App import successful")',
        ],
        appWorkingDir
      );
    }

    console.log('\n=== Test 6: Gunicorn Startup Simulation ===');
    await testGunicornHealth(pythonExecutablePath, packagedPythonSourceDir);

    console.log(
      '\n✓ All tests passed! The Python executable should work correctly.'
    );
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
      env: {
        ...createPythonEnv(),
      },
    });

    let stdout = '';
    let stderr = '';
    let resolved = false;

    // Set timeout
    const timer = setTimeout(() => {
      if (!resolved) {
        childProcess.kill('SIGTERM');
        resolved = true;
        reject(new Error(`Command timed out after ${timeout}ms`));
      }
    }, timeout);

    childProcess.stdout?.on('data', data => {
      const output = data.toString();
      stdout += output;
      console.log(`STDOUT: ${output.trim()}`);

      // For Gunicorn, if we see it starting to listen, that's success
      if (
        output.includes('Listening at:') ||
        output.includes('Booting worker')
      ) {
        if (!resolved) {
          console.log('✓ Gunicorn startup detected - killing process');
          clearTimeout(timer);
          childProcess.kill('SIGTERM');
          resolved = true;
          resolve();
        }
      }
    });

    childProcess.stderr?.on('data', data => {
      const output = data.toString();
      stderr += output;
      console.log(`STDERR: ${output.trim()}`);
    });

    childProcess.on('close', code => {
      clearTimeout(timer);
      if (!resolved) {
        if (code === 0) {
          console.log('✓ Command completed successfully');
          resolve();
        } else {
          console.log(`✗ Command failed with exit code: ${code}`);
          if (stderr) console.log(`Error output: ${stderr}`);
          reject(
            new Error(
              `Command failed with exit code ${code}: ${stderr || 'No error message'}`
            )
          );
        }
        resolved = true;
      }
    });

    childProcess.on('error', error => {
      clearTimeout(timer);
      if (!resolved) {
        console.log(`✗ Process error: ${error.message}`);
        reject(error);
        resolved = true;
      }
    });
  });
}

function getFreePort() {
  return new Promise((resolve, reject) => {
    const server = net.createServer();
    server.listen(0, '127.0.0.1', () => {
      const address = server.address();
      server.close(() => {
        if (address && typeof address === 'object') {
          resolve(address.port);
        } else {
          reject(new Error('Failed to allocate a local port'));
        }
      });
    });
    server.on('error', reject);
  });
}

function waitForHealth(port, timeout = 30000) {
  const deadline = Date.now() + timeout;
  const url = `http://127.0.0.1:${port}/health`;

  return new Promise((resolve, reject) => {
    const check = () => {
      const request = http.get(
        url,
        { headers: { 'X-Auth-Token': AUTH_TOKEN } },
        response => {
          response.resume();
          if (response.statusCode === 200) {
            resolve();
            return;
          }
          reject(
            new Error(`Health check returned HTTP ${response.statusCode}`)
          );
        }
      );

      request.on('error', error => {
        if (Date.now() >= deadline) {
          reject(new Error(`Health check timed out: ${error.message}`));
          return;
        }
        setTimeout(check, 500);
      });
      request.setTimeout(2000, () => {
        request.destroy(new Error('Health check request timed out'));
      });
    };

    check();
  });
}

function waitForClose(childProcess, timeout = 5000) {
  return new Promise(resolve => {
    let resolved = false;
    const timer = setTimeout(() => {
      if (!resolved) {
        childProcess.kill('SIGKILL');
      }
    }, timeout);

    childProcess.on('close', (code, signal) => {
      if (!resolved) {
        resolved = true;
        clearTimeout(timer);
        resolve({ code, signal });
      }
    });
  });
}

async function testGunicornHealth(pythonPath, workDir) {
  const port = await getFreePort();
  const gunicornArgs = [
    '-m',
    'gunicorn',
    '--bind',
    `127.0.0.1:${port}`,
    '--workers',
    '1',
    '--worker-class',
    'uvicorn.workers.UvicornWorker',
    '--timeout',
    '30',
    '--log-level',
    'info',
    'app:app',
  ];

  console.log(`Command: ${pythonPath} ${gunicornArgs.join(' ')}`);
  console.log(`Working directory: ${workDir}`);
  console.log(`Health check URL: http://127.0.0.1:${port}/health`);

  const childProcess = spawn(pythonPath, gunicornArgs, {
    cwd: workDir,
    stdio: ['pipe', 'pipe', 'pipe'],
    env: createPythonEnv(port),
  });

  childProcess.stdout?.on('data', data => {
    console.log(`STDOUT: ${data.toString().trim()}`);
  });
  childProcess.stderr?.on('data', data => {
    console.log(`STDERR: ${data.toString().trim()}`);
  });

  const closePromise = new Promise(resolve => {
    childProcess.on('close', (code, signal) => {
      resolve({ code, signal });
    });
  });

  try {
    await Promise.race([
      waitForHealth(port),
      closePromise.then(({ code, signal }) => {
        throw new Error(
          `Gunicorn exited before health check passed (code=${code}, signal=${signal || 'none'})`
        );
      }),
    ]);
    console.log('✓ Gunicorn health check passed with X-Auth-Token');
  } finally {
    if (childProcess.exitCode === null) {
      const closeAfterKill = waitForClose(childProcess);
      childProcess.kill('SIGTERM');
      await closeAfterKill;
    }
  }
}
