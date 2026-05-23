#!/usr/bin/env node

/**
 * Docker run script for building with existing image
 * Handles cross-platform path resolution for Docker volumes
 */

const { execSync } = require('child_process');
const path = require('path');
const fs = require('fs');

// Get absolute path to project root
const projectRoot = process.cwd();
const distPath = path.join(projectRoot, 'dist');
const releasePath = path.join(projectRoot, 'release');
const condaEnvPath = path.join(projectRoot, 'conda_env');
const pythonDistPath = path.join(projectRoot, 'python_dist');

// Convert Windows path to Docker-compatible format if needed
function dockerizePath(filePath) {
  if (process.platform === 'win32') {
    // Convert backslashes to forward slashes for Docker
    return filePath.replace(/\\/g, '/');
  }
  return filePath;
}

const dockerDistPath = dockerizePath(distPath);
const dockerReleasePath = dockerizePath(releasePath);
const dockerCondaEnvPath = dockerizePath(condaEnvPath);
const dockerPythonDistPath = dockerizePath(pythonDistPath);

const args = process.argv.slice(2);
const useDryRun = args.includes('--dry-run');

console.log('Running Docker build with existing image...');
console.log(`Dist path: ${dockerDistPath}`);
console.log(`Release path: ${dockerReleasePath}`);
console.log(`Conda env path: ${dockerCondaEnvPath}`);
console.log(`Python dist path: ${dockerPythonDistPath}`);

const dockerCmd = `docker run --rm -v "${dockerDistPath}:/app/dist" -v "${dockerReleasePath}:/app/release" -v "${dockerCondaEnvPath}:/app/conda_env" -v "${dockerPythonDistPath}:/app/python_dist" pyscf-front-builder`;

if (useDryRun) {
  console.log('\n=== Dry run: command was not executed ===');
  console.log(`Run command: ${dockerCmd}`);
  process.exit(0);
}

try {
  fs.mkdirSync(distPath, { recursive: true });
  fs.mkdirSync(releasePath, { recursive: true });

  // Run container with volume mounts
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
