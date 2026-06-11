# Security Best Practices Report

## Executive Summary

No critical or high-severity application vulnerabilities were identified in this source review. The project already has several strong controls: a random loopback backend auth token, strict CORS/TrustedHost settings, request size limits, redacted server errors, parameterized SQLite queries, Electron navigation guards, IPC sender validation, disabled Node integration, context isolation, CSP, and restrictive permissions for local secret-bearing files.

The main residual risks are Electron-specific hardening gaps and local-secret handling: untrusted chat-rendered links can reach `shell.openExternal`, packaged windows still use `file://`, CSP allows all loopback ports, the Gemini API key is stored as plaintext JSON protected only by file permissions, and runtime GPU package installation retains supply-chain risk. These are most relevant if renderer JavaScript is compromised through a future XSS, vulnerable dependency, or malicious AI-generated content.

## Scope

Reviewed stack:

- Backend: Python 3.12, FastAPI, Gunicorn/Uvicorn, SQLite, requests.
- Frontend: TypeScript, React 19, React Markdown, Electron renderer.
- Desktop shell: Electron 42, preload IPC bridge, main-process backend orchestration.

Assumptions:

- This is a local desktop app with a loopback FastAPI backend, not a public multi-user web service.
- The backend auth token is intended to protect local HTTP endpoints from unrelated local pages/processes.
- Findings are based on source review and targeted grep/static checks, not dynamic exploit testing.

Guidance consulted:

- Local skill references:
  - `python-fastapi-web-server-security.md`
  - `javascript-general-web-frontend-security.md`
  - `javascript-typescript-react-web-frontend-security.md`
- Electron official security guidance:
  - https://www.electronjs.org/docs/latest/tutorial/security
  - https://www.electronjs.org/docs/latest/api/structures/web-preferences

## Positive Controls Observed

- FastAPI docs are disabled on the main app (`src/python/app.py:390-391`), while development docs are only registered outside packaged mode (`src/python/api/__init__.py:30-35`) and use pinned Swagger UI SRI hashes plus CSP (`src/python/api/swagger_ui.py:20-57`, `src/python/api/swagger_ui.py:83-108`).
- Backend auth requires `X-Auth-Token` except narrowly scoped development docs/OPTIONS/test bypasses (`src/python/app.py:169-207`).
- Electron generates a random 32-byte auth token and passes it only to the backend/main-process network injector (`src/main.ts:50-53`, `src/main/python-env.ts:272-289`, `src/main/backend-auth-injection.ts:61-85`).
- CORS is restricted to loopback dev origins or packaged `file://`/`null` origins with `allow_credentials=False` (`src/python/app.py:210-246`).
- TrustedHostMiddleware is configured with loopback-only hosts in non-dev environments (`src/python/app.py:406-430`).
- Request body size enforcement rejects oversized or missing `Content-Length` body-bearing requests (`src/python/app.py:96-166`).
- Unknown 5xx details are redacted before returning to clients (`src/python/app.py:296-355`).
- Electron windows set `nodeIntegration: false` and `contextIsolation: true` (`src/main/window-manager.ts:33-38`, `src/main/splash-window-manager.ts:32-37`).
- Electron permission requests are denied by default, CSP is injected for HTTP content, and navigation/window-open guards are installed (`src/main/session-hardening.ts:77-118`, `src/main/renderer-entry.ts:220-247`).
- Privileged IPC handlers validate the sender (`src/main/ipc.ts:13-84`, `src/main/ipc-security.ts:13-29`).
- SQLite queries use parameterized placeholders for user-controlled values (`src/python/database/chat_history.py:313-319`, `src/python/database/chat_history.py:342-349`, `src/python/database/chat_history.py:448-450`, `src/python/database/chat_history.py:488-503`).
- Chat history DB and settings files are hardened with owner-only POSIX permissions (`src/python/database/chat_history.py:19-22`, `src/python/database/chat_history.py:140-169`, `src/python/quantum_calc/settings_manager.py:15-18`, `src/python/quantum_calc/settings_manager.py:244-264`).

## Critical Findings

None identified.

## High Findings

None identified.

## Medium Findings

### M-001: Untrusted chat-rendered URLs can reach `shell.openExternal`

- Rule ID: ELECTRON-OPENEXTERNAL-001
- Severity: Medium
- Location:
  - `src/web/components/ChatMessage.tsx:26-36`
  - `src/preload.ts:21-28`
  - `src/main/ipc-security.ts:35-87`
  - `src/main/ipc.ts:33-45`
- Evidence:
  - Model-generated markdown links are rendered and HTTP(S) links call `window.electronAPI.openExternalUrl(href)` from the renderer.
  - The preload exposes `openExternalUrl`.
  - Main-process validation limits URLs to non-empty HTTP(S), length <= 2048, and a parseable URL, then calls `shell.openExternal`.
- Impact:
  - A malicious AI response or any future renderer compromise can cause arbitrary HTTP(S) URLs to be opened in the user's default browser. Protocol filtering significantly reduces host command-execution risk, but this still enables phishing and limited data exfiltration via URL parameters.
- Fix:
  - Treat model-generated links as untrusted content. Prefer an interstitial confirmation showing the full destination, or enforce an allowlist for known trusted domains.
  - Add explicit blocking for localhost/private IP/link-local destinations unless intentionally needed.
  - Consider passing a user-gesture token from the click handler to the IPC call so arbitrary renderer code cannot invoke `openExternalUrl` silently.
- Mitigation:
  - Current HTTP(S)-only validation and IPC sender checks reduce the risk; keep them.
- False positive notes:
  - This is not currently an arbitrary-protocol issue because `file:`, custom protocols, and non-HTTP(S) schemes are rejected.

### M-002: Packaged renderer uses `file://` instead of a custom app protocol

- Rule ID: ELECTRON-FILE-PROTOCOL-001
- Severity: Medium
- Location:
  - `src/main/renderer-entry.ts:120-132`
  - `src/main/renderer-entry.ts:164-170`
  - `src/main/window-manager.ts:109-114`
  - `src/main/splash-window-manager.ts:53-57`
- Evidence:
  - Packaged entries are returned as `type: 'file'`.
  - Main and splash windows call `loadFile(...)`.
  - Navigation guards allow `file:` URLs within the renderer directory.
- Impact:
  - Electron's own guidance recommends custom protocols over `file://` because `file://` has special local-file behavior. If renderer JavaScript is ever compromised, `file://` increases the blast radius compared with an app-scoped custom protocol.
- Fix:
  - Register a privileged custom protocol such as `app://` with `protocol.handle`.
  - Serve only packaged renderer assets through that protocol.
  - Move CSP delivery to response headers for the custom protocol and remove broad `file://` navigation allowance.
- Mitigation:
  - Existing navigation guards restrict same-window navigation to the renderer directory, and CSP is present in HTML/meta plus HTTP headers where applicable.
- False positive notes:
  - No current XSS path was found in this review; this is defense-in-depth for future renderer compromise.

### M-003: CSP allows renderer connections to every loopback port

- Rule ID: ELECTRON-CSP-CONNECT-001
- Severity: Medium
- Location:
  - `index.html:10-12`
  - `src/main/session-hardening.ts:35-47`
  - `src/main/session-hardening.ts:56-68`
- Evidence:
  - `connect-src` allows `http://127.0.0.1:*` and `ws://127.0.0.1:*`.
- Impact:
  - If renderer script execution is compromised, the page can attempt requests to arbitrary local services on 127.0.0.1. CORS may prevent reading many responses, but state-changing requests or services with permissive CORS remain exposed.
- Fix:
  - Generate packaged CSP with the actual backend port instead of a wildcard.
  - If using a custom protocol, deliver CSP as a dynamic response header from the main process.
  - Keep the wildcard only in development if needed for Vite/HMR and document it as dev-only.
- Mitigation:
  - Backend requests still require `X-Auth-Token`, and the token injector only targets the selected backend port (`src/main/backend-auth-injection.ts:61-85`).
- False positive notes:
  - The wildcard is useful for a dynamic backend port, so the fix should preserve runtime port flexibility rather than hard-coding port 5000.

### M-004: Gemini API key is stored as plaintext JSON at rest

- Rule ID: SECRET-STORAGE-001
- Severity: Medium
- Location:
  - `src/python/quantum_calc/settings_manager.py:123-133`
  - `src/python/quantum_calc/settings_manager.py:226-267`
  - `src/python/services/settings_service.py:24-32`
  - `src/python/api/settings.py:19-36`
- Evidence:
  - `gemini_api_key` is part of persisted app settings.
  - `save_settings` writes the full settings dict to JSON.
  - API responses mask the key, but internal settings keep the plaintext value for Gemini calls.
- Impact:
  - Any local process or user account that can read the settings file can recover the Gemini API key. POSIX `0600` permissions reduce accidental exposure but do not provide credential-store protection.
- Fix:
  - Store the Gemini key in the OS credential store/keychain and keep only a configured flag or key reference in JSON.
  - If cross-platform keychain integration is deferred, document the current local-storage threat model explicitly in settings UI/docs.
- Mitigation:
  - Current HTTP responses mask the key and logs use `mask_settings`; POSIX permissions are hardened to `0700`/`0600`.
- False positive notes:
  - This is a local desktop app, so plaintext-at-rest may be an accepted tradeoff during development.

### M-005: Runtime GPU package installation keeps residual supply-chain risk

- Rule ID: DEP-SUPPLYCHAIN-001
- Severity: Medium
- Location:
  - `src/python/services/system_service.py:224-255`
  - `src/python/services/system_service.py:450-488`
  - `src/python/api/system.py:111-140`
- Evidence:
  - The app can run `python -m pip install ...` at runtime for GPU4PySCF packages.
  - The code explicitly notes that transitive dependencies are not hash-pinned.
  - Production builds deny runtime install unless `PYSCF_ALLOW_RUNTIME_INSTALL=1`; development/test allow it.
- Impact:
  - Runtime package installation expands the trusted supply chain after packaging and may install changed transitive artifacts if upstream package indexes change or are compromised.
- Fix:
  - Prefer prebuilt, locked GPU environments per supported CUDA generation.
  - If runtime install must remain, generate a hash-locked wheelhouse/requirements set for each CUDA family and install from that controlled source.
  - Keep `PYSCF_ALLOW_RUNTIME_INSTALL` fail-closed in production.
- Mitigation:
  - The endpoint requires normal backend auth, local loopback source, explicit confirmation, Linux support, pinned top-level package versions, and no `shell=True`.
- False positive notes:
  - This is a deliberate product feature, not accidental command injection.

## Low Findings

### L-001: PubChem path segments are built from raw query text

- Rule ID: OUTBOUND-URL-001
- Severity: Low
- Location:
  - `src/python/api/pubchem.py:16-30`
  - `src/python/services/pubchem_service.py:17-24`
  - `src/python/pubchem/client.py:80-89`
- Evidence:
  - `search_type` is allowlisted, but `query.strip()` is interpolated directly into the PubChem URL path.
- Impact:
  - The host is fixed to PubChem, so this is not SSRF into arbitrary infrastructure. However, special characters such as `/`, `?`, and `#` can alter the intended PubChem path/query semantics and cause unexpected upstream requests.
- Fix:
  - Percent-encode user query path segments, for example with `urllib.parse.quote(query.strip(), safe="")`.
  - Add tests for spaces, slashes, `?`, `#`, and non-ASCII compound names.
- Mitigation:
  - Query length is capped at the API boundary (`MAX_QUERY_LENGTH = 500`), and errors avoid reflecting proprietary query strings.
- False positive notes:
  - Requests may normalize some characters internally, but explicit path encoding is safer and clearer.

### L-002: Dependency vulnerability scanning was not executable in this environment

- Rule ID: DEP-SCANNING-001
- Severity: Low
- Location:
  - `package.json:20-45`
  - `.github/environment.yml:12-37`
- Evidence:
  - `npm audit --audit-level=moderate --omit=dev --json` failed because the sandbox could not resolve `registry.npmjs.org`.
  - `~/miniforge3/envs/pyscf-env/bin/python -m pip_audit --format json` failed because `pip_audit` is not installed.
- Impact:
  - The source review cannot confirm whether current npm or Python dependencies have known vulnerabilities.
- Fix:
  - Add a CI job for `npm audit` or an equivalent lockfile scanner.
  - Add `pip-audit` or an equivalent Python/conda vulnerability scanner to CI, ideally using the generated lock files.
- Mitigation:
  - Dependencies are mostly version-pinned or lockfile-backed (`package-lock.json`, `.github/pyscf-env.conda-lock.yml`, explicit conda lock outputs).
- False positive notes:
  - This is a verification gap, not proof of a vulnerable dependency.

## Verification Performed

- Loaded project rules from `.claude/rules/*` and `.claude/docs/DESIGN.md`.
- Loaded security guidance for FastAPI, general frontend JavaScript/TypeScript, and React.
- Inspected Electron/FastAPI/auth/CORS/CSP/IPC/file/path/SQL/secret/storage code with targeted `rg` and line-numbered reads.
- Checked git status before writing this report; the working tree was clean.
- Ran `npm audit --audit-level=moderate --omit=dev --json`; it failed due DNS/network restriction to `registry.npmjs.org`.
- Ran `~/miniforge3/envs/pyscf-env/bin/python -m pip_audit --format json`; it failed because `pip_audit` is not installed.

## Recommended Fix Order

1. Address `M-001` first because it directly exposes an OS-level capability to renderer-triggered untrusted URLs.
2. Address `M-003` together with `M-002`; a custom protocol makes dynamic, stricter CSP header delivery much easier.
3. Move `gemini_api_key` to OS credential storage when the settings subsystem is next touched.
4. Decide whether runtime GPU installation remains a product requirement; if yes, document the risk and move toward locked wheels.
5. Percent-encode PubChem path segments as a small, low-risk cleanup.
