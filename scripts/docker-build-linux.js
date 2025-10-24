#!/usr/bin/env node

/**
 * Docker build script for Linux distribution
 * Handles cross-platform path resolution for Docker volumes
 */

const { execSync } = require('child_process');
const path = require('path');
const fs = require('fs');

// Get absolute path to project root
const projectRoot = process.cwd();
const distPath = path.join(projectRoot, 'dist');

// Convert Windows path to Docker-compatible format if needed
function dockerizePath(filePath) {
  if (process.platform === 'win32') {
    // Convert C:\Users\... to /c/Users/... or use Windows path directly
    // Docker Desktop on Windows handles Windows paths
    return filePath.replace(/\\/g, '/');
  }
  return filePath;
}

const dockerDistPath = dockerizePath(distPath);

console.log('Building Docker image for Linux package...');
console.log(`Project root: ${projectRoot}`);
console.log(`Dist path: ${dockerDistPath}`);

try {
  // Clean dist directory if it exists
  if (fs.existsSync(distPath)) {
    console.log('\n=== Step 0: Cleaning dist directory ===');
    try {
      // Try to remove with retries (helps with OneDrive sync issues)
      fs.rmSync(distPath, { recursive: true, force: true, maxRetries: 3, retryDelay: 100 });
      console.log('Cleaned dist directory');
    } catch (cleanError) {
      console.warn('⚠️  Warning: Could not fully clean dist directory:');
      console.warn(cleanError.message);
      console.warn('\nThis is usually caused by OneDrive sync or file locks.');
      console.warn('Solutions:');
      console.warn('  1. Manually delete the dist folder in File Explorer');
      console.warn('  2. Pause OneDrive sync temporarily');
      console.warn('  3. Move the project outside OneDrive');
      console.warn('\nPress Ctrl+C to cancel, or we will continue in 5 seconds...\n');

      // Give user time to cancel if they want to try manual cleanup
      execSync('timeout /t 5 /nobreak', { stdio: 'inherit' });
    }
  }

  // Build Docker image
  console.log('\n=== Step 1: Building Docker image ===');
  execSync('docker build -t pyscf-front-builder .', {
    stdio: 'inherit',
    cwd: projectRoot
  });

  console.log('\n=== Step 2: Running build in container ===');
  console.log(`Mounting: ${dockerDistPath} -> /app/dist`);

  // Run container with volume mount
  const dockerCmd = `docker run --rm -v "${dockerDistPath}:/app/dist" pyscf-front-builder`;

  execSync(dockerCmd, {
    stdio: 'inherit',
    cwd: projectRoot
  });

  console.log('\n✅ Linux package built successfully!');
  console.log(`Output: ${distPath}`);
} catch (error) {
  console.error('\n❌ Build failed:', error.message);
  process.exit(1);
}
