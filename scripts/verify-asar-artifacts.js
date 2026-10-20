#!/usr/bin/env node

const fs = require('fs');
const path = require('path');
const asar = require('@electron/asar');

const DEFAULT_ROOTS = ['release', 'dist'];
const MAX_SAMPLE = 80;

const forbiddenArtifacts = [
  { name: 'DMG artifact', pattern: /(^|\/)[^/]+\.dmg(?:\.blockmap)?$/i },
  {
    name: 'AppImage artifact',
    pattern: /(^|\/)[^/]+\.AppImage(?:\.blockmap)?$/i,
  },
  { name: 'ZIP artifact', pattern: /(^|\/)[^/]+\.zip(?:\.blockmap)?$/i },
  { name: 'mac unpacked directory', pattern: /(^|\/)mac-[^/]+(?=\/|$)/i },
  {
    name: 'linux unpacked directory',
    pattern: /(^|\/)linux-unpacked(?=\/|$)/i,
  },
  {
    name: 'unpacked package directory',
    pattern: /(^|\/)[^/]+-unpacked(?=\/|$)/i,
  },
  {
    name: 'electron-builder debug file',
    pattern: /(^|\/)builder-debug\.yml$/i,
  },
];

function findAsarFiles(directory) {
  if (!fs.existsSync(directory)) {
    return [];
  }

  const entries = fs.readdirSync(directory, { withFileTypes: true });
  const archives = [];

  for (const entry of entries) {
    const entryPath = path.join(directory, entry.name);

    if (entry.isDirectory()) {
      archives.push(...findAsarFiles(entryPath));
      continue;
    }

    if (entry.isFile() && entry.name === 'app.asar') {
      archives.push(entryPath);
    }
  }

  return archives;
}

function findForbiddenEntries(archivePath) {
  return asar
    .listPackage(archivePath)
    .map(entry => entry.replace(/^\/+/, ''))
    .filter(entry =>
      forbiddenArtifacts.some(({ pattern }) => pattern.test(entry))
    );
}

const roots = process.argv.slice(2);
const scanRoots = roots.length > 0 ? roots : DEFAULT_ROOTS;
const archives = scanRoots.flatMap(findAsarFiles);

if (archives.length === 0) {
  console.log(
    'No app.asar files found under ' +
      scanRoots.join(', ') +
      '; skipping artifact contamination check.'
  );
  process.exit(0);
}

let hasForbiddenEntries = false;

for (const archivePath of archives) {
  const forbiddenEntries = findForbiddenEntries(archivePath);

  if (forbiddenEntries.length === 0) {
    console.log(
      'OK: ' + archivePath + ' does not contain packaging artifacts.'
    );
    continue;
  }

  hasForbiddenEntries = true;
  console.error(
    'ERROR: ' +
      archivePath +
      ' contains ' +
      forbiddenEntries.length +
      ' packaging artifact entries.'
  );

  for (const entry of forbiddenEntries.slice(0, MAX_SAMPLE)) {
    console.error('  - ' + entry);
  }

  if (forbiddenEntries.length > MAX_SAMPLE) {
    console.error(
      '  ... ' +
        (forbiddenEntries.length - MAX_SAMPLE) +
        ' more entries omitted'
    );
  }
}

if (hasForbiddenEntries) {
  process.exit(1);
}

console.log(
  'Checked ' +
    archives.length +
    ' app.asar file(s); no packaging artifacts found.'
);
