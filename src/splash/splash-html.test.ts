const assert: typeof import('node:assert/strict') = require('node:assert/strict');
const fs: typeof import('node:fs') = require('node:fs');
const path: typeof import('node:path') = require('node:path');

const splashHtmlPath = path.resolve(__dirname, '../../splash.html');

const readSplashHtml = (): string => {
  return fs.readFileSync(splashHtmlPath, 'utf8');
};

const testSplashHtmlIncludesCriticalStylesBeforeModuleScript = (): void => {
  const html = readSplashHtml();
  const criticalStyleIndex = html.indexOf('data-splash-critical');
  const moduleScriptIndex = html.indexOf('<script type="module"');

  assert.notEqual(criticalStyleIndex, -1);
  assert.notEqual(moduleScriptIndex, -1);
  assert.ok(criticalStyleIndex < moduleScriptIndex);
};

const run = (): void => {
  testSplashHtmlIncludesCriticalStylesBeforeModuleScript();
  console.log('splash html tests passed');
};

run();

export {};
