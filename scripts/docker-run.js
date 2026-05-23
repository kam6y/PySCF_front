#!/usr/bin/env node

/**
 * Docker run script for building with an existing image.
 * Builds inside /app, then copies artifacts to host-mounted output folders.
 */

const { execSync } = require('child_process');
const path = require('path');
const fs = require('fs');

const IMAGE_NAME = 'pyscf-front-builder';
const DOCKER_PLATFORM = 'linux/amd64';
const CONTAINER_OUTPUT_ROOT = '/host-output';

const projectRoot = process.cwd();
const distPath = path.join(projectRoot, 'dist');
const releasePath = path.join(projectRoot, 'release');
const condaEnvPath = path.join(projectRoot, 'conda_env');
const pythonDistPath = path.join(projectRoot, 'python_dist');

function dockerizePath(filePath) {
  if (process.platform === 'win32') {
    return filePath.replace(/\\/g, '/');
  }
  return filePath;
}

function shellQuote(value) {
  if (/^[A-Za-z0-9_/:=.,@%+-]+$/.test(value)) {
    return value;
  }

  const escaped = value.replace(/'/g, "'\\''");
  return `'${escaped}'`;
}

function commandText(args) {
  return args.map(shellQuote).join(' ');
}

function ensureOutputDirectory(directoryPath) {
  fs.mkdirSync(directoryPath, { recursive: true });
}


const dockerDistPath = dockerizePath(distPath);
const dockerReleasePath = dockerizePath(releasePath);
const dockerCondaEnvPath = dockerizePath(condaEnvPath);
const dockerPythonDistPath = dockerizePath(pythonDistPath);

const args = process.argv.slice(2);
const useDryRun = args.includes('--dry-run');

console.log('Running Docker build with existing image...');
console.log(`Docker platform: ${DOCKER_PLATFORM}`);
console.log(`Dist output path: ${dockerDistPath}`);
console.log(`Release output path: ${dockerReleasePath}`);
console.log(`Conda env output path: ${dockerCondaEnvPath}`);
console.log(`Python dist output path: ${dockerPythonDistPath}`);

const containerBuildScript = [
  'source /root/miniforge3/etc/profile.d/conda.sh',
  'conda activate pyscf-env',
  'npm run package:linux',
  `mkdir -p ${CONTAINER_OUTPUT_ROOT}/dist ${CONTAINER_OUTPUT_ROOT}/release ${CONTAINER_OUTPUT_ROOT}/conda_env ${CONTAINER_OUTPUT_ROOT}/python_dist`,
  `find ${CONTAINER_OUTPUT_ROOT}/dist -mindepth 1 -maxdepth 1 -exec rm -rf {} +`,
  `find ${CONTAINER_OUTPUT_ROOT}/release -mindepth 1 -maxdepth 1 -exec rm -rf {} +`,
  `find ${CONTAINER_OUTPUT_ROOT}/conda_env -mindepth 1 -maxdepth 1 -exec rm -rf {} +`,
  `find ${CONTAINER_OUTPUT_ROOT}/python_dist -mindepth 1 -maxdepth 1 -exec rm -rf {} +`,
  `cp -a /app/dist/. ${CONTAINER_OUTPUT_ROOT}/dist/`,
  `cp -a /app/release/. ${CONTAINER_OUTPUT_ROOT}/release/`,
  `cp -a /app/conda_env/. ${CONTAINER_OUTPUT_ROOT}/conda_env/`,
  `if [ -d /app/python_dist ]; then cp -a /app/python_dist/. ${CONTAINER_OUTPUT_ROOT}/python_dist/; fi`,
].join(' && ');

const dockerCmd = commandText([
  'docker',
  'run',
  '--platform',
  DOCKER_PLATFORM,
  '--rm',
  '-v',
  `${dockerDistPath}:${CONTAINER_OUTPUT_ROOT}/dist`,
  '-v',
  `${dockerReleasePath}:${CONTAINER_OUTPUT_ROOT}/release`,
  '-v',
  `${dockerCondaEnvPath}:${CONTAINER_OUTPUT_ROOT}/conda_env`,
  '-v',
  `${dockerPythonDistPath}:${CONTAINER_OUTPUT_ROOT}/python_dist`,
  IMAGE_NAME,
  'bash',
  '-lc',
  containerBuildScript,
]);

if (useDryRun) {
  console.log('\n=== Dry run: command was not executed ===');
  console.log(`Run command: ${dockerCmd}`);
  console.log('\nContainer output copy:');
  console.log(`  clean ${CONTAINER_OUTPUT_ROOT}/dist/*, then /app/dist/. -> ${CONTAINER_OUTPUT_ROOT}/dist/`);
  console.log(`  clean ${CONTAINER_OUTPUT_ROOT}/release/*, then /app/release/. -> ${CONTAINER_OUTPUT_ROOT}/release/`);
  console.log(`  clean ${CONTAINER_OUTPUT_ROOT}/conda_env/*, then /app/conda_env/. -> ${CONTAINER_OUTPUT_ROOT}/conda_env/`);
  console.log(`  clean ${CONTAINER_OUTPUT_ROOT}/python_dist/*, then /app/python_dist/. -> ${CONTAINER_OUTPUT_ROOT}/python_dist/ (if present)`);
  process.exit(0);
}

try {
  ensureOutputDirectory(distPath);
  ensureOutputDirectory(releasePath);
  ensureOutputDirectory(condaEnvPath);
  ensureOutputDirectory(pythonDistPath);

  console.log('\n=== Running build in container ===');
  console.log(`Mounting: ${dockerDistPath} -> ${CONTAINER_OUTPUT_ROOT}/dist`);
  console.log(`Mounting: ${dockerReleasePath} -> ${CONTAINER_OUTPUT_ROOT}/release`);
  console.log(`Mounting: ${dockerCondaEnvPath} -> ${CONTAINER_OUTPUT_ROOT}/conda_env`);
  console.log(`Mounting: ${dockerPythonDistPath} -> ${CONTAINER_OUTPUT_ROOT}/python_dist`);
  console.log('Container cleanup runs only against internal /app paths.');

  execSync(dockerCmd, {
    stdio: 'inherit',
    cwd: projectRoot,
  });

  console.log('\n✅ Build completed successfully!');
  console.log(`Output: ${releasePath}`);
} catch (error) {
  console.error('\n❌ Build failed:', error.message);
  process.exit(1);
}
