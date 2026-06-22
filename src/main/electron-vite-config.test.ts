const assert: typeof import('node:assert/strict') = require('node:assert/strict');
const fs: typeof import('node:fs') = require('node:fs');
const path: typeof import('node:path') = require('node:path');
const ts: typeof import('typescript') = require('typescript');

const configPath = path.resolve(__dirname, '../../electron.vite.config.ts');

const getObjectProperty = (
  objectLiteral: import('typescript').ObjectLiteralExpression,
  propertyName: string
): import('typescript').Expression | null => {
  const property = objectLiteral.properties.find(
    (
      item
    ): item is import('typescript').PropertyAssignment & {
      name: import('typescript').Identifier;
    } => {
      return (
        ts.isPropertyAssignment(item) &&
        ts.isIdentifier(item.name) &&
        item.name.text === propertyName
      );
    }
  );

  return property?.initializer ?? null;
};

const readRendererOptimizeDepsInclude = (): string[] => {
  const sourceText = fs.readFileSync(configPath, 'utf8');
  const sourceFile = ts.createSourceFile(
    configPath,
    sourceText,
    ts.ScriptTarget.Latest,
    true
  );

  const rendererConfigDeclaration = sourceFile.statements.find(statement => {
    return (
      ts.isVariableStatement(statement) &&
      statement.declarationList.declarations.some(
        declaration =>
          ts.isIdentifier(declaration.name) &&
          declaration.name.text === 'rendererConfig'
      )
    );
  });

  assert.ok(rendererConfigDeclaration);
  assert.ok(ts.isVariableStatement(rendererConfigDeclaration));

  const rendererConfig =
    rendererConfigDeclaration.declarationList.declarations.find(
      declaration =>
        ts.isIdentifier(declaration.name) &&
        declaration.name.text === 'rendererConfig'
    )?.initializer;

  assert.ok(rendererConfig);
  assert.ok(ts.isObjectLiteralExpression(rendererConfig));

  const optimizeDeps = getObjectProperty(rendererConfig, 'optimizeDeps');
  assert.ok(optimizeDeps);
  assert.ok(ts.isObjectLiteralExpression(optimizeDeps));

  const include = getObjectProperty(optimizeDeps, 'include');
  assert.ok(include);
  assert.ok(ts.isArrayLiteralExpression(include));

  return include.elements.map(element => {
    assert.ok(ts.isStringLiteral(element));
    return element.text;
  });
};

const readRendererOptimizeDepsExclude = (): string[] => {
  const sourceText = fs.readFileSync(configPath, 'utf8');
  const sourceFile = ts.createSourceFile(
    configPath,
    sourceText,
    ts.ScriptTarget.Latest,
    true
  );

  const rendererConfigDeclaration = sourceFile.statements.find(statement => {
    return (
      ts.isVariableStatement(statement) &&
      statement.declarationList.declarations.some(
        declaration =>
          ts.isIdentifier(declaration.name) &&
          declaration.name.text === 'rendererConfig'
      )
    );
  });

  assert.ok(rendererConfigDeclaration);
  assert.ok(ts.isVariableStatement(rendererConfigDeclaration));

  const rendererConfig =
    rendererConfigDeclaration.declarationList.declarations.find(
      declaration =>
        ts.isIdentifier(declaration.name) &&
        declaration.name.text === 'rendererConfig'
    )?.initializer;

  assert.ok(rendererConfig);
  assert.ok(ts.isObjectLiteralExpression(rendererConfig));

  const optimizeDeps = getObjectProperty(rendererConfig, 'optimizeDeps');
  assert.ok(optimizeDeps);
  assert.ok(ts.isObjectLiteralExpression(optimizeDeps));

  const exclude = getObjectProperty(optimizeDeps, 'exclude');
  if (!exclude) {
    return [];
  }

  assert.ok(ts.isArrayLiteralExpression(exclude));

  return exclude.elements.map(element => {
    assert.ok(ts.isStringLiteral(element));
    return element.text;
  });
};

const testKetcherCommonJsDependenciesArePrebundled = (): void => {
  const include = readRendererOptimizeDepsInclude();

  assert.ok(include.includes('hoist-non-react-statics'));
  assert.ok(include.includes('ketcher-core'));
  assert.ok(include.includes('ketcher-react'));
  assert.ok(include.includes('ketcher-standalone'));
  assert.ok(include.includes('lodash'));
  assert.ok(include.includes('lodash/fp'));
};

const testKetcherReactIsNotExcludedFromPrebundling = (): void => {
  const exclude = readRendererOptimizeDepsExclude();

  assert.ok(!exclude.includes('ketcher-react'));
};

// ============================================================
// J1/J13: devCspRelaxPlugin tests — simulate transformIndexHtml
// ============================================================

/**
 * Simulate the Vite CSP plugins against the actual HTML files to detect marker drift.
 */

const indexHtmlPath = path.resolve(__dirname, '../../index.html');
const splashHtmlPath = path.resolve(__dirname, '../../splash.html');

const readHtml = (htmlPath: string): string => fs.readFileSync(htmlPath, 'utf8');

/**
 * Simulate the devCspRelaxPlugin transformIndexHtml logic.
 * We replicate the exact string replacements from the plugin source
 * to verify they work on the actual HTML content.
 */
const simulateDevCspRelax = (html: string): string => {
  let result = html.replace(
    "script-src 'self' 'wasm-unsafe-eval'",
    "script-src 'self' 'wasm-unsafe-eval' 'unsafe-eval' 'unsafe-inline'"
  );
  if (result === html) {
    throw new Error('script-src marker not found');
  }

  const beforeConnect = result;
  result = result.replace(
    "connect-src 'self'",
    "connect-src 'self' http://127.0.0.1:* ws://127.0.0.1:*"
  );
  if (result === beforeConnect) {
    throw new Error('connect-src marker not found');
  }

  return result;
};

/**
 * Simulate the stripMetaCspPlugin transformIndexHtml logic.
 */
const simulateStripMetaCsp = (html: string): string => {
  const metaCspPattern = /\s*<meta\s[\s\S]*?http-equiv=["']Content-Security-Policy["'][\s\S]*?\/?>\s*/i;
  const result = html.replace(metaCspPattern, '\n');
  if (result === html) {
    throw new Error('meta CSP tag not found');
  }
  return result;
};

// J1: devCspRelaxPlugin works on index.html
const testDevCspRelax_indexHtml = (): void => {
  const html = readHtml(indexHtmlPath);
  const relaxed = simulateDevCspRelax(html);

  // Must have relaxed script-src
  assert.ok(
    relaxed.includes("'unsafe-eval'"),
    'index.html: relaxed CSP must include unsafe-eval'
  );
  assert.ok(
    relaxed.includes("'unsafe-inline'"),
    'index.html: relaxed CSP must include unsafe-inline in script-src'
  );
  // Must have relaxed connect-src
  assert.ok(
    relaxed.includes('http://127.0.0.1:*'),
    'index.html: relaxed CSP must include connect-src wildcard'
  );
};

// J1: devCspRelaxPlugin works on splash.html (the regression case)
const testDevCspRelax_splashHtml = (): void => {
  const html = readHtml(splashHtmlPath);
  // This would throw before the J1 fix (splash.html lacked wasm-unsafe-eval)
  const relaxed = simulateDevCspRelax(html);

  assert.ok(
    relaxed.includes("'unsafe-eval'"),
    'splash.html: relaxed CSP must include unsafe-eval'
  );
  assert.ok(
    relaxed.includes('http://127.0.0.1:*'),
    'splash.html: relaxed CSP must include connect-src wildcard'
  );
};

// J13: devCspRelaxPlugin throws when marker is genuinely absent
const testDevCspRelax_failLoudOnMissingMarker = (): void => {
  assert.throws(
    () => simulateDevCspRelax('<html><head></head><body></body></html>'),
    /script-src marker not found/,
    'Must throw when CSP marker is absent'
  );
};

// J13: stripMetaCspPlugin strips meta CSP from index.html
const testStripMetaCsp_indexHtml = (): void => {
  const html = readHtml(indexHtmlPath);
  const stripped = simulateStripMetaCsp(html);

  assert.ok(
    !stripped.includes('Content-Security-Policy'),
    'index.html: meta CSP tag must be stripped'
  );
  // Other content must survive
  assert.ok(
    stripped.includes('<div id="root">'),
    'index.html: non-CSP content must survive'
  );
};

// J13: stripMetaCspPlugin strips meta CSP from splash.html
const testStripMetaCsp_splashHtml = (): void => {
  const html = readHtml(splashHtmlPath);
  const stripped = simulateStripMetaCsp(html);

  assert.ok(
    !stripped.includes('Content-Security-Policy'),
    'splash.html: meta CSP tag must be stripped'
  );
  // Other content must survive
  assert.ok(
    stripped.includes('splash-container'),
    'splash.html: non-CSP content must survive'
  );
};

// J13: stripMetaCspPlugin throws when tag is absent
const testStripMetaCsp_failLoudOnMissingTag = (): void => {
  assert.throws(
    () => simulateStripMetaCsp('<html><head></head><body></body></html>'),
    /meta CSP tag not found/,
    'Must throw when meta CSP tag is absent'
  );
};

const run = (): void => {
  testKetcherCommonJsDependenciesArePrebundled();
  testKetcherReactIsNotExcludedFromPrebundling();

  // J1/J13: devCspRelaxPlugin
  testDevCspRelax_indexHtml();
  testDevCspRelax_splashHtml();
  testDevCspRelax_failLoudOnMissingMarker();

  // J13: stripMetaCspPlugin
  testStripMetaCsp_indexHtml();
  testStripMetaCsp_splashHtml();
  testStripMetaCsp_failLoudOnMissingTag();

  console.log('electron-vite config tests passed (8 tests)');
};

run();

export {};
