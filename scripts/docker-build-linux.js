#!/usr/bin/env node

/**
 * Docker build script for Linux distribution.
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

const args = process.argv.slice(2);
const useNoCache = args.includes('--no-cache');
const useClean = args.includes('--clean');

console.log('Building Docker image for Linux package...');
console.log(`Project root: ${projectRoot}`);
console.log(`Docker platform: ${DOCKER_PLATFORM}`);
console.log(`Dist output path: ${dockerDistPath}`);
console.log(`Release output path: ${dockerReleasePath}`);

if (process.platform === 'win32' && projectRoot.includes('OneDrive')) {
  console.warn('\n⚠️  WARNING: Project is located in OneDrive folder');
  console.warn('This may cause build issues due to file sync conflicts.');
  console.warn(
    'Consider moving the project outside OneDrive for better stability.\n'
  );
}

const buildFlags = [];
if (useNoCache) {
  buildFlags.push('--no-cache');
}

const containerBuildScript = [
  'source /root/miniforge3/etc/profile.d/conda.sh',
  'conda activate pyscf-env',
  'npm run package:linux',
  `mkdir -p ${CONTAINER_OUTPUT_ROOT}/dist ${CONTAINER_OUTPUT_ROOT}/release`,
  `find ${CONTAINER_OUTPUT_ROOT}/dist -mindepth 1 -maxdepth 1 -exec rm -rf {} +`,
  `find ${CONTAINER_OUTPUT_ROOT}/release -mindepth 1 -maxdepth 1 -exec rm -rf {} +`,
  `cp -a /app/dist/. ${CONTAINER_OUTPUT_ROOT}/dist/`,
  `cp -a /app/release/. ${CONTAINER_OUTPUT_ROOT}/release/`,
].join(' && ');

const buildCmd = commandText([
  'docker',
  'build',
  '--platform',
  DOCKER_PLATFORM,
  ...buildFlags,
  '-t',
  IMAGE_NAME,
  '.',
]);

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
  IMAGE_NAME,
  'bash',
  '-lc',
  containerBuildScript,
]);

try {
  if (useClean) {
    console.log('\n=== Step 0a: Cleaning Docker build cache ===');
    try {
      execSync('docker builder prune -f', { stdio: 'inherit' });
      console.log('Docker build cache cleaned');
    } catch (error) {
      console.warn('⚠️  Warning: Could not clean Docker cache:', error.message);
    }
  }

  ensureOutputDirectory(distPath);
  ensureOutputDirectory(releasePath);

  console.log('\n=== Step 1: Building Docker image ===');
  if (useNoCache) {
    console.log('Using --no-cache flag (slower but avoids cache issues)');
  }
  console.log(`Command: ${buildCmd}\n`);

  execSync(buildCmd, {
    stdio: 'inherit',
    cwd: projectRoot,
    env: {
      ...process.env,
      DOCKER_BUILDKIT: '0',
    },
  });

  console.log('\n=== Step 2: Running build in container ===');
  console.log(`Mounting: ${dockerDistPath} -> ${CONTAINER_OUTPUT_ROOT}/dist`);
  console.log(
    `Mounting: ${dockerReleasePath} -> ${CONTAINER_OUTPUT_ROOT}/release`
  );
  console.log('Container cleanup runs only against internal /app paths.');

  execSync(dockerCmd, {
    stdio: 'inherit',
    cwd: projectRoot,
  });

  console.log('\n✅ Linux package built successfully!');
  console.log(`Output: ${releasePath}`);
} catch (error) {
  console.error('\n❌ Build failed:', error.message);
  process.exit(1);
}
