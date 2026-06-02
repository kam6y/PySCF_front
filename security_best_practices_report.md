# Security Best Practices Report

## Executive Summary

This audit reviewed the Electron + React + FastAPI/Python application against the project rules, the security-best-practices skill guidance for FastAPI and React, and high-signal Electron hardening checks. The strongest issue is local secret persistence: the Gemini API key is stored and returned in plaintext, and observed local settings/chat files are world-readable under the current umask. No critical remote unauthenticated issue was confirmed; backend auth is generated with cryptographic randomness, SQL usage in chat history is parameterized, and calculation ID path resolution rejects traversal.

Recommended priority:

1. Protect local secrets and chat data at rest, and stop returning the full API key from `GET /api/settings`.
2. Add body/schema limits for large text fields and expensive endpoints.
3. Harden the packaged renderer CSP and production Electron surface.
4. Rework the runtime GPU4PySCF installer into a pinned/confirmed/offline-safe flow.
5. Self-host or pin development Swagger UI assets.

## Scope and Method

Reviewed:

- FastAPI entrypoint, middleware, routers, settings, chat, system/GPU install, and calculation APIs.
- React renderer markdown/link handling, runtime API clients, CSP in HTML entrypoints, local storage use.
- Electron main/preload IPC boundaries, BrowserWindow webPreferences, menu exposure, backend token generation.
- Dependency and CI posture from `package.json`, lockfiles, `.github/environment.yml`, and GitHub workflows.

Commands used included `rg` pattern scans, `nl -ba` line inspection, `stat` permission checks, and official Electron documentation lookup for current renderer sandbox/security behavior.

## Positive Observations

- `src/main.ts:50-53` generates a per-run 32-byte random backend auth token with `crypto.randomBytes`.
- `src/main/window-manager.ts:28-32` disables Node integration and enables context isolation for the main renderer.
- `src/python/app.py:202` disables FastAPI's default docs/OpenAPI routes, and `src/python/api/__init__.py:30-35` keeps custom Swagger UI out of packaged mode.
- `src/python/quantum_calc/_calculation_repository.py:34-53` rejects path separators and verifies calculation paths stay inside the base directory.
- `src/python/database/chat_history.py:193-196`, `src/python/database/chat_history.py:219-223`, and related DB calls use SQLite parameters.

## High Severity

### SEC-001: Local secrets and chat data are persisted with default world-readable permissions

Rule ID: GENERAL-SECRETS-001, REACT-AUTH-001, FASTAPI-RESP-001

Severity: High

Location:

- `src/python/quantum_calc/settings_manager.py:40-45`
- `src/python/quantum_calc/settings_manager.py:191-195`
- `src/python/services/settings_service.py:111-127`
- `src/python/api/settings.py:11-14`
- `src/web/pages/SettingsPage.tsx:84-100`
- `src/python/database/chat_history.py:81-96`

Evidence:

```python
# src/python/quantum_calc/settings_manager.py:40-45
app_data_dir = base_dir / ".pyscf_native_app"
app_data_dir.mkdir(exist_ok=True)
self.settings_file = app_data_dir / "app-settings.json"
```

```python
# src/python/quantum_calc/settings_manager.py:191-195
with open(temp_file, 'w', encoding='utf-8') as f:
    json.dump(settings_dict, f, indent=2, ensure_ascii=False)
temp_file.replace(self.settings_file)
```

```python
# src/python/services/settings_service.py:124-127
settings = get_current_settings()
return settings.model_dump(mode='json')
```

The audit also observed these permissions without reading file contents:

```text
-rw-r--r-- /Users/goodapple/.pyscf_native_app/app-settings.json
drwxr-xr-x /Users/goodapple/.pyscf_native_app
-rw-r--r-- /Users/goodapple/.pyscf_app/data/chat_history.db
drwxr-xr-x /Users/goodapple/.pyscf_app/data
```

Impact:

If a Gemini API key is configured, another local OS user or backup/indexing process with filesystem access can read it from the settings JSON; chat history can also contain sensitive prompts/results and is similarly readable.

Fix:

- Create secret/data directories with `0o700`, and create/replace secret-bearing files with `0o600`.
- Migrate existing files/directories by `chmod` on startup.
- Prefer OS keychain storage for `gemini_api_key` instead of JSON.
- Change `GET /api/settings` to return `gemini_api_key_configured: boolean` or a masked value, and make the API key write-only except when intentionally replacing it.
- Disable production DevTools or at least avoid putting plaintext secrets in renderer state.

Mitigation:

Warn users not to store a Gemini key until the storage migration lands; manually run `chmod 700 ~/.pyscf_native_app ~/.pyscf_app/data` and `chmod 600 ~/.pyscf_native_app/app-settings.json ~/.pyscf_app/data/chat_history.db` on existing installations.

False positive notes:

This is most severe on multi-user machines or managed workstations. Same-user malware can usually read user-owned files regardless of mode bits, so keychain storage is still preferable.

## Medium Severity

### SEC-002: Large request bodies and expensive user-controlled inputs lack consistent size caps

Rule ID: FASTAPI-LIMITS-001, FASTAPI-VALID-001, REACT-NET-001

Severity: Medium

Location:

- `src/python/api/quantum.py:30-39`
- `src/python/api/system.py:77-99`
- `src/python/api/agent.py:20-22`
- `src/python/api/agent.py:96-105`
- `src/api-spec/openapi.yaml:1396-1399`
- `src/api-spec/openapi.yaml:1445-1448`
- `src/api-spec/openapi.yaml:3795-3821`

Evidence:

```python
# src/python/api/quantum.py:30-39
raw_body = await request.body()
if not raw_body:
    raise HTTPException(status_code=400, detail='Request body is required')

raw_data = json.loads(raw_body)
```

```yaml
# src/api-spec/openapi.yaml:1396-1399
xyz:
  type: string
  minLength: 1
  description: XYZ molecular structure data
```

```yaml
# src/api-spec/openapi.yaml:3799-3817
message:
  type: string
history:
  type: array
  items:
    type: object
```

Impact:

An authenticated caller, compromised renderer, or development-mode local webpage can send very large JSON, XYZ, Ketcher data, or chat history payloads that consume memory, trigger expensive quantum work, or increase LLM API cost.

Fix:

- Add a central FastAPI request-size middleware that rejects oversized bodies with `413`.
- Add OpenAPI/Pydantic `maxLength` for `xyz`, `ketcher_data`, SMILES/PubChem query text, chat `message`, history item text, and session names where missing.
- Add `maxItems` for chat history and molecular atom arrays where appropriate.
- Lower `MAX_MESSAGE_LENGTH` from `100000` unless that size is a deliberate product requirement, and cap total history characters before sending to Gemini.
- Add focused tests for oversized quantum, agent, and GPU install payloads.

Mitigation:

Keep the backend bound to `127.0.0.1` and require `PYSCF_AUTH_TOKEN` outside development, which the current packaged app already does.

False positive notes:

Some large molecules legitimately need bigger XYZ data. The fix should use chemistry-aware limits rather than a tiny generic cap.

### SEC-003: Packaged renderer CSP still allows `unsafe-eval`

Rule ID: REACT-CSP-001, JS-CSP-002

Severity: Medium

Location:

- `index.html:5-9`

Evidence:

```html
<meta
  http-equiv="Content-Security-Policy"
  content="default-src 'self'; script-src 'self' 'unsafe-eval'; style-src 'self' 'unsafe-inline'; connect-src 'self' http://127.0.0.1:* ws://127.0.0.1:*; img-src 'self' data:; worker-src 'self' blob:; child-src 'self' blob:;"
/>
```

Impact:

If an XSS or dependency-rendering bug is introduced, `unsafe-eval` weakens CSP's ability to block string-to-code execution in the Electron renderer.

Fix:

- Remove `unsafe-eval` from the production `index.html` CSP.
- If Vite/development tooling requires `unsafe-eval`, generate separate dev and packaged CSPs.
- Consider Trusted Types in report-only mode for renderer DOM sinks after the CSP is tightened.

Mitigation:

The current scan did not find direct `eval`, `new Function`, `dangerouslySetInnerHTML`, or `innerHTML` usage in renderer code, so the primary issue is defense-in-depth rather than a confirmed exploit path.

False positive notes:

Some 3D/chemistry libraries may require eval-like behavior. If so, document the exact dependency and constrain the policy to development or an isolated viewer context.

### SEC-004: Runtime GPU4PySCF installation uses pip against live package indexes with unpinned fallback packages

Rule ID: FASTAPI-INJECT-002, FASTAPI-SUPPLY-001

Severity: Medium

Location:

- `src/python/api/system.py:77-104`
- `src/python/services/system_service.py:148-167`
- `src/python/services/system_service.py:191-215`
- `src/python/services/system_service.py:325-380`

Evidence:

```python
# src/python/services/system_service.py:154-159
pip_command = [sys.executable, "-m", "pip", "install", "--no-cache-dir", "--prefer-binary"]
...
pip_command.extend(packages)
```

```python
# src/python/services/system_service.py:209-215
deps_latest = [
    f"cupy-cuda{cuda_suffix}",
    libxc_requirement,
]
if include_cutensor:
    deps_latest.append(f"cutensor-cu{cuda_major}")
```

Impact:

When invoked, the application mutates its Python environment by downloading/installing executable Python packages at runtime; unpinned fallback packages expand supply-chain risk and can change behavior across installs.

Fix:

- Prefer packaging GPU4PySCF variants through the existing locked build/release process.
- If runtime install remains required, pin exact versions, use hashes or a trusted internal wheelhouse, and show an explicit user confirmation before running installation.
- Return only sanitized install diagnostics; avoid exposing long pip output in UI/logs if it can contain local paths or environment details.

Mitigation:

The endpoint requires the app auth token and also checks that the client address is loopback in `src/python/api/system.py:80-92`, and package names are selected from fixed maps rather than raw request strings.

False positive notes:

This is a deliberate feature, not shell injection: `subprocess.run` uses an argument list and `shell=True` was not found.

## Low Severity

### SEC-005: Development Swagger UI loads third-party assets without SRI or a strict CSP

Rule ID: REACT-SRI-001, JS-SUPPLY-001, FASTAPI-OPENAPI-001

Severity: Low

Location:

- `src/python/app.py:62-70`
- `src/python/api/__init__.py:30-35`
- `src/python/api/swagger_ui.py:10-29`

Evidence:

```python
# src/python/app.py:62-70
def _is_development_api_docs_path(path: str) -> bool:
    return path == '/api-docs' or path.startswith('/api-docs/')

if _is_development_api_docs_path(request.url.path) and os.getenv('PYSCF_ENV') == 'development':
    return await call_next(request)
```

```html
<!-- src/python/api/swagger_ui.py:18-24 -->
<link rel="stylesheet" href="https://unpkg.com/swagger-ui-dist/swagger-ui.css">
<script src="https://unpkg.com/swagger-ui-dist/swagger-ui-bundle.js"></script>
<script>
  SwaggerUIBundle({ url: '/api-docs/spec.json', dom_id: '#swagger-ui' });
</script>
```

Impact:

The development docs page runs remote JavaScript from `unpkg.com` with no integrity pinning, and it is intentionally reachable without the app token in development.

Fix:

- Self-host `swagger-ui-dist` as a pinned dependency, or add exact-version URLs with SRI.
- Add a CSP for the docs response and remove inline script by moving initialization to a static JS file.
- Consider making unauthenticated dev docs opt-in through a dedicated environment flag.

Mitigation:

The router is not registered in packaged mode (`src/python/api/__init__.py:30-35`), so the risk is development-only.

False positive notes:

This is not a production exposure unless `PYSCF_RESOURCES_PATH` packaging detection is bypassed or the development server is intentionally exposed.

### SEC-006: Production Electron window leaves DevTools available

Rule ID: ELECTRON-HARDENING-001

Severity: Low

Location:

- `src/main/window-manager.ts:28-33`
- `src/main/menu.ts:108-115`

Evidence:

```ts
// src/main/window-manager.ts:28-33
webPreferences: {
  preload: path.join(__dirname, 'preload.js'),
  nodeIntegration: false,
  contextIsolation: true,
  devTools: true,
},
```

```ts
// src/main/menu.ts:108-115
{
  label: 'View',
  submenu: [
    { role: 'reload' as const },
    { role: 'forceReload' as const },
    { role: 'toggleDevTools' as const },
```

Impact:

Any local user can inspect renderer state and call exposed preload APIs more easily; this compounds SEC-001 because the current renderer receives the Gemini API key.

Fix:

- Set `devTools: !app.isPackaged` for the main window.
- Remove `toggleDevTools`, `reload`, and `forceReload` menu items in packaged builds unless needed for support diagnostics.
- If support diagnostics require DevTools, gate it behind an explicit debug flag.

Mitigation:

Node integration is disabled and context isolation is enabled, which are the most important Electron renderer safeguards.

False positive notes:

For a local scientific desktop application, production DevTools may be an intentional support tradeoff. Document that decision if retained.

## Verification Notes

- No direct renderer use of `dangerouslySetInnerHTML`, `.innerHTML`, `eval`, `new Function`, broad `postMessage`, or token-bearing `localStorage` was found in the scanned source.
- No `shell=True` subprocess calls were found in source paths reviewed.
- SQL query construction in chat history uses parameters.
- `npm ci` and conda lockfiles are used in CI/release workflows; no online vulnerability audit was run in this pass.
- FastAPI/Starlette versions in the lockfiles are newer than the historical Starlette examples referenced by the skill guidance.

## Recommended Fix Order

1. Fix SEC-001 first: file modes, keychain or masked settings API, and migration for existing files.
2. Add request-size middleware and schema limits for SEC-002.
3. Split dev/packaged CSP and remove production `unsafe-eval`.
4. Decide whether runtime GPU package installation remains a product requirement; if yes, pin/hash it.
5. Self-host Swagger UI assets and restrict development docs behavior.
6. Disable production DevTools unless there is an explicit support requirement.

## Sources

- Electron Security Checklist: https://www.electronjs.org/docs/latest/tutorial/security
- Electron WebPreferences, including sandbox defaults: https://www.electronjs.org/docs/latest/api/structures/web-preferences
- Electron Process Sandboxing: https://www.electronjs.org/docs/latest/tutorial/sandbox/
- Local skill references used:
  - `python-fastapi-web-server-security.md`
  - `javascript-typescript-react-web-frontend-security.md`
  - `javascript-general-web-frontend-security.md`
