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

const run = (): void => {
  testKetcherCommonJsDependenciesArePrebundled();
  testKetcherReactIsNotExcludedFromPrebundling();
  console.log('electron-vite config tests passed');
};

run();

export {};
