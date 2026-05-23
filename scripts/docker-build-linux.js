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
const releasePath = path.join(projectRoot, 'release');

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
const dockerReleasePath = dockerizePath(releasePath);

// Parse command line arguments
const args = process.argv.slice(2);
const useNoCache = args.includes('--no-cache');
const useClean = args.includes('--clean');
const useDryRun = args.includes('--dry-run');

console.log('Building Docker image for Linux package...');
console.log(`Project root: ${projectRoot}`);
console.log(`Dist path: ${dockerDistPath}`);
console.log(`Release path: ${dockerReleasePath}`);

// Check if project is in OneDrive (Windows only)
if (process.platform === 'win32' && projectRoot.includes('OneDrive')) {
  console.warn('\n⚠️  WARNING: Project is located in OneDrive folder');
  console.warn('This may cause build issues due to file sync conflicts.');
  console.warn(
    'Consider moving the project outside OneDrive for better stability.\n'
  );
}

// Construct commands with options
const buildFlags = [];
if (useNoCache) {
  buildFlags.push('--no-cache');
}

const buildFlagText = buildFlags.length > 0 ? `${buildFlags.join(' ')} ` : '';
const buildCmd = `docker build ${buildFlagText}-t pyscf-front-builder .`;
const dockerCmd = `docker run --rm -v "${dockerDistPath}:/app/dist" -v "${dockerReleasePath}:/app/release" pyscf-front-builder`;

if (useDryRun) {
  console.log('\n=== Dry run: commands were not executed ===');
  console.log(`Build command: ${buildCmd}`);
  console.log(`Run command: ${dockerCmd}`);
  process.exit(0);
}

try {
  // Clean Docker build cache if requested
  if (useClean) {
    console.log('\n=== Step 0a: Cleaning Docker build cache ===');
    try {
      execSync('docker builder prune -f', { stdio: 'inherit' });
      console.log('Docker build cache cleaned');
    } catch (error) {
      console.warn('⚠️  Warning: Could not clean Docker cache:', error.message);
    }
  }

  // Clean dist directory if it exists
  if (fs.existsSync(distPath)) {
    console.log('\n=== Step 0b: Cleaning dist directory ===');
    try {
      // Try to remove with retries (helps with OneDrive sync issues)
      fs.rmSync(distPath, {
        recursive: true,
        force: true,
        maxRetries: 3,
        retryDelay: 100,
      });
      console.log('Cleaned dist directory');
    } catch (cleanError) {
      console.warn('⚠️  Warning: Could not fully clean dist directory:');
      console.warn(cleanError.message);
      console.warn('\nThis is usually caused by OneDrive sync or file locks.');
      console.warn('Solutions:');
      console.warn('  1. Manually delete the dist folder in File Explorer');
      console.warn('  2. Pause OneDrive sync temporarily');
      console.warn('  3. Move the project outside OneDrive');
      console.warn('  4. Run with --clean flag to clear Docker cache');
      console.warn(
        '\nPress Ctrl+C to cancel, or we will continue in 5 seconds...\n'
      );

      // Give user time to cancel if they want to try manual cleanup
      execSync('timeout /t 5 /nobreak', { stdio: 'inherit' });
    }
  }

  fs.mkdirSync(distPath, { recursive: true });
  fs.mkdirSync(releasePath, { recursive: true });

  // Build Docker image
  console.log('\n=== Step 1: Building Docker image ===');

  if (useNoCache) {
    console.log('Using --no-cache flag (slower but avoids cache issues)');
  }
  console.log(`Command: ${buildCmd}\n`);

  // Use legacy builder (DOCKER_BUILDKIT=0) to avoid BuildKit cache issues
  execSync(buildCmd, {
    stdio: 'inherit',
    cwd: projectRoot,
    env: {
      ...process.env,
      DOCKER_BUILDKIT: '0', // Disable BuildKit to avoid cache corruption issues
    },
  });

  console.log('\n=== Step 2: Running build in container ===');
  console.log(`Mounting: ${dockerDistPath} -> /app/dist`);
  console.log(`Mounting: ${dockerReleasePath} -> /app/release`);

  // Run container with volume mounts
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
