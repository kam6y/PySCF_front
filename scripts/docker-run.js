#!/usr/bin/env node

/**
 * Docker run script for building with existing image
 * Handles cross-platform path resolution for Docker volumes
 */

const { execSync } = require('child_process');
const path = require('path');

// Get absolute path to project root
const projectRoot = process.cwd();
const distPath = path.join(projectRoot, 'dist');
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
const dockerCondaEnvPath = dockerizePath(condaEnvPath);
const dockerPythonDistPath = dockerizePath(pythonDistPath);

console.log('Running Docker build with existing image...');
console.log(`Dist path: ${dockerDistPath}`);
console.log(`Conda env path: ${dockerCondaEnvPath}`);
console.log(`Python dist path: ${dockerPythonDistPath}`);

try {
  // Run container with volume mounts
  const dockerCmd = `docker run --rm -v "${dockerDistPath}:/app/dist" -v "${dockerCondaEnvPath}:/app/conda_env" -v "${dockerPythonDistPath}:/app/python_dist" pyscf-front-builder`;

  execSync(dockerCmd, {
    stdio: 'inherit',
    cwd: projectRoot
  });

  console.log('\n✅ Build completed successfully!');
  console.log(`Output: ${distPath}`);
} catch (error) {
  console.error('\n❌ Build failed:', error.message);
  process.exit(1);
}
