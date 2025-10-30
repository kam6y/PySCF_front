#!/usr/bin/env node

/**
 * Version bump script for PySCF Native App
 *
 * This script updates version across all project files and creates a git tag.
 *
 * Usage:
 *   node scripts/bump-version.js <new-version> [--no-push] [--force]
 *   npm run release -- <new-version> [--no-push] [--force]
 *
 * Examples:
 *   npm run release -- 0.0.2
 *   npm run release -- 0.1.0 --no-push
 *   npm run release -- 0.0.1 --force  # Force recreate existing tag
 */

const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');

// ANSI color codes for terminal output
const colors = {
  reset: '\x1b[0m',
  green: '\x1b[32m',
  yellow: '\x1b[33m',
  red: '\x1b[31m',
  blue: '\x1b[34m',
  cyan: '\x1b[36m'
};

function log(message, color = 'reset') {
  console.log(`${colors[color]}${message}${colors.reset}`);
}

function error(message) {
  log(`❌ ${message}`, 'red');
}

function success(message) {
  log(`✅ ${message}`, 'green');
}

function info(message) {
  log(`ℹ️  ${message}`, 'cyan');
}

function warning(message) {
  log(`⚠️  ${message}`, 'yellow');
}

// Parse command line arguments
const args = process.argv.slice(2);
const newVersion = args.find(arg => !arg.startsWith('--'));
const noPush = args.includes('--no-push');
const forceTag = args.includes('--force');

if (!newVersion) {
  error('Usage: npm run release -- <new-version> [--no-push] [--force]');
  error('Example: npm run release -- 0.0.2');
  error('Options:');
  error('  --no-push  Create commit and tag locally without pushing');
  error('  --force    Force recreate tag if it already exists');
  process.exit(1);
}

// Validate version format (semver-like)
const versionRegex = /^\d+\.\d+\.\d+(-[a-zA-Z0-9.-]+)?$/;
if (!versionRegex.test(newVersion)) {
  error(`Invalid version format: ${newVersion}`);
  error('Version must follow semver format: X.Y.Z or X.Y.Z-label');
  error('Examples: 0.0.2, 0.1.0, 1.0.0, 0.0.2-beta.1');
  process.exit(1);
}

// Get project root directory
const projectRoot = path.resolve(__dirname, '..');

// Define files to update
const filesToUpdate = [
  {
    path: 'package.json',
    update: (content) => {
      const pkg = JSON.parse(content);
      const oldVersion = pkg.version;
      pkg.version = newVersion;
      info(`package.json: ${oldVersion} → ${newVersion}`);
      return JSON.stringify(pkg, null, 2) + '\n';
    }
  },
  {
    path: 'src/api-spec/openapi.yaml',
    update: (content) => {
      const versionLine = content.match(/version: [\d.]+(-[a-zA-Z0-9.-]+)?/);
      if (versionLine) {
        const oldVersion = versionLine[0].split(': ')[1];
        info(`openapi.yaml: ${oldVersion} → ${newVersion}`);
      }
      return content.replace(
        /version: [\d.]+(-[a-zA-Z0-9.-]+)?/,
        `version: ${newVersion}`
      );
    }
  },
  {
    path: 'src/python/__init__.py',
    update: (content) => {
      const versionLine = content.match(/__version__ = ["'][\d.]+(-[a-zA-Z0-9.-]+)?["']/);
      if (versionLine) {
        const oldVersion = versionLine[0].match(/["'](.*?)["']/)[1];
        info(`__init__.py: ${oldVersion} → ${newVersion}`);
      }
      return content.replace(
        /__version__ = ["'][\d.]+(-[a-zA-Z0-9.-]+)?["']/,
        `__version__ = "${newVersion}"`
      );
    }
  },
  {
    path: 'config/server-config.json',
    update: (content) => {
      const config = JSON.parse(content);
      const oldVersion = config.app_info?.version;
      if (oldVersion) {
        info(`server-config.json: ${oldVersion} → ${newVersion}`);
      }
      config.app_info = config.app_info || {};
      config.app_info.version = newVersion;
      return JSON.stringify(config, null, 2) + '\n';
    }
  },
  {
    path: 'src/python/config.py',
    update: (content) => {
      const versionMatch = content.match(/"version": "[\d.]+(-[a-zA-Z0-9.-]+)?"/);
      if (versionMatch) {
        const oldVersion = versionMatch[0].match(/"([\d.]+(-[a-zA-Z0-9.-]+)?)"/)[1];
        info(`config.py (fallback): ${oldVersion} → ${newVersion}`);
      }
      return content.replace(
        /"version": "[\d.]+(-[a-zA-Z0-9.-]+)?"/,
        `"version": "${newVersion}"`
      );
    }
  }
];

function updateVersionFiles() {
  log('\n📝 Updating version files...', 'blue');

  filesToUpdate.forEach(({ path: filePath, update }) => {
    const fullPath = path.join(projectRoot, filePath);

    if (!fs.existsSync(fullPath)) {
      warning(`File not found: ${filePath}`);
      return;
    }

    try {
      const content = fs.readFileSync(fullPath, 'utf8');
      const updatedContent = update(content);
      fs.writeFileSync(fullPath, updatedContent, 'utf8');
    } catch (err) {
      error(`Failed to update ${filePath}: ${err.message}`);
      process.exit(1);
    }
  });

  success('All version files updated!');
}

function checkGitStatus() {
  log('\n🔍 Checking git status...', 'blue');

  try {
    const status = execSync('git status --porcelain', { encoding: 'utf8' });
    const uncommittedChanges = status
      .split('\n')
      .filter(line => line.trim())
      .filter(line => {
        // Ignore changes to version files
        return !filesToUpdate.some(f => line.includes(f.path));
      });

    if (uncommittedChanges.length > 0) {
      warning('You have uncommitted changes:');
      uncommittedChanges.forEach(line => console.log(`  ${line}`));
      warning('Please commit or stash these changes before creating a release.');
      process.exit(1);
    }

    success('Git working directory is clean');
  } catch (err) {
    error(`Failed to check git status: ${err.message}`);
    process.exit(1);
  }
}

function getCurrentBranch() {
  try {
    return execSync('git rev-parse --abbrev-ref HEAD', { encoding: 'utf8' }).trim();
  } catch (err) {
    error(`Failed to get current branch: ${err.message}`);
    process.exit(1);
  }
}

function checkExistingTag() {
  log('\n🏷️  Checking for existing tag...', 'blue');

  const tagName = `v${newVersion}`;
  let localTagExists = false;
  let remoteTagExists = false;

  // Check local tag
  try {
    execSync(`git rev-parse ${tagName}`, { encoding: 'utf8', stdio: 'pipe' });
    localTagExists = true;
  } catch (err) {
    // Tag doesn't exist locally
  }

  // Check remote tag
  try {
    execSync(`git ls-remote --tags origin ${tagName}`, { encoding: 'utf8', stdio: 'pipe' });
    const output = execSync(`git ls-remote --tags origin ${tagName}`, { encoding: 'utf8' });
    if (output.trim()) {
      remoteTagExists = true;
    }
  } catch (err) {
    // Remote check failed, probably no remote
  }

  if (localTagExists || remoteTagExists) {
    warning(`Tag ${tagName} already exists!`);
    if (localTagExists) {
      warning(`  - Found in local repository`);
    }
    if (remoteTagExists) {
      warning(`  - Found in remote repository`);
    }

    if (!forceTag) {
      error('');
      error('To recreate this tag, use the --force option:');
      error(`  npm run release -- ${newVersion} --force`);
      process.exit(1);
    }

    log('\n🔨 Force mode enabled - will recreate tag', 'yellow');

    // Delete local tag if it exists
    if (localTagExists) {
      try {
        execSync(`git tag -d ${tagName}`, { stdio: 'inherit' });
        success(`Deleted local tag ${tagName}`);
      } catch (err) {
        error(`Failed to delete local tag: ${err.message}`);
        process.exit(1);
      }
    }

    // Delete remote tag if it exists
    if (remoteTagExists && !noPush) {
      try {
        execSync(`git push origin :refs/tags/${tagName}`, { stdio: 'inherit' });
        success(`Deleted remote tag ${tagName}`);
      } catch (err) {
        error(`Failed to delete remote tag: ${err.message}`);
        process.exit(1);
      }
    }
  } else {
    success(`Tag ${tagName} does not exist - ready to create`);
  }
}

function commitAndTag() {
  log('\n📦 Creating git commit and tag...', 'blue');

  const currentBranch = getCurrentBranch();
  info(`Current branch: ${currentBranch}`);

  try {
    // Stage version files
    const filePaths = filesToUpdate.map(f => f.path).join(' ');
    execSync(`git add ${filePaths}`, { stdio: 'inherit' });

    // Create commit
    const commitMessage = `Release v${newVersion}

- すべてのバージョンを${newVersion}に更新

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>`;

    execSync(`git commit -m "${commitMessage.replace(/"/g, '\\"')}"`, { stdio: 'inherit' });
    success(`Commit created: Release v${newVersion}`);

    // Create tag
    execSync(`git tag v${newVersion}`, { stdio: 'inherit' });
    success(`Tag created: v${newVersion}`);

    // Push if not disabled
    if (!noPush) {
      log('\n🚀 Pushing to remote...', 'blue');
      execSync(`git push origin ${currentBranch}`, { stdio: 'inherit' });
      execSync(`git push origin v${newVersion}`, { stdio: 'inherit' });
      success('Pushed to remote repository');

      log('\n🎉 Release completed!', 'green');
      info(`GitHub Actions will now build and create release v${newVersion}`);
      info(`Check progress at: https://github.com/kam6y/PySCF_front/actions`);
      info(`Release will be available at: https://github.com/kam6y/PySCF_front/releases/tag/v${newVersion}`);
    } else {
      log('\n✋ Skipping push (--no-push flag)', 'yellow');
      warning('Remember to push manually:');
      info(`  git push origin ${currentBranch}`);
      info(`  git push origin v${newVersion}`);
    }
  } catch (err) {
    error(`Failed to commit and tag: ${err.message}`);
    process.exit(1);
  }
}

// Main execution
function main() {
  log('\n🚀 PySCF Native App - Version Bump Script', 'cyan');
  log('==========================================\n', 'cyan');

  info(`New version: ${newVersion}`);
  info(`Push to remote: ${!noPush ? 'Yes' : 'No (--no-push flag)'}`);
  info(`Force mode: ${forceTag ? 'Yes (--force flag)' : 'No'}\n`);

  checkGitStatus();
  checkExistingTag();
  updateVersionFiles();
  commitAndTag();
}

main();
