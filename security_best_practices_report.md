# Security Best Practices Report

作成日: 2026-06-09

## Executive Summary

このリポジトリは Electron + React/TypeScript + FastAPI/Python のローカルデスクトップアプリとして、認証トークン、CORS 制限、TrustedHostMiddleware、CSP、Electron のナビゲーションガード、IPC 送信元検証など、主要な防御策がすでに入っています。

監査で Critical は見つかりませんでした。優先度が最も高いのは、利用中の Electron 37 系が公式サポート終了済みで、Chromium/Node.js のセキュリティ修正を取り逃がすリスクです。次点で、CI の依存関係スキャン不足と、ロック済み Python 環境を CI 中に直接 `pip install` で変更している点が供給網リスクとして残ります。

## Scope And Method

対象:

- Electron main/preload/splash: `src/main*`, `src/preload.ts`, `src/splash/*`
- React/TypeScript renderer: `src/web/*`, `index.html`, `splash.html`
- FastAPI/Python backend: `src/python/*`
- API contract and generated models: `src/api-spec/openapi.yaml`, `src/python/generated_models.py`, `src/web/types/generated-api.ts`
- CI and dependency metadata: `package.json`, `package-lock.json`, `.github/workflows/*`, `.github/environment.yml`, `.github/*conda-lock*`

参照したガイダンス:

- `security-best-practices` skill references for FastAPI, React, and frontend JavaScript/TypeScript.
- Electron official security checklist and release/support policy.
- react-markdown official documentation via Context7 for raw HTML and URL transform behavior.

実行した主な確認:

- `rg` による危険 sink 検索: `dangerouslySetInnerHTML`, `innerHTML`, `eval`, `new Function`, `postMessage`, `localStorage`, `shell.openExternal`, `FileResponse`, `subprocess`, `CORSMiddleware`, `TrustedHostMiddleware` など。
- `package-lock.json` 上の Electron 実解決バージョン確認。
- `npm audit --audit-level=high --omit=dev --json` はネットワーク制限により失敗。

## Positive Security Controls Confirmed

- Electron renderer は `nodeIntegration: false` / `contextIsolation: true` で作成されています: `src/main/window-manager.ts:33`, `src/main/window-manager.ts:35`, `src/main/window-manager.ts:36`, `src/main/splash-window-manager.ts:32`, `src/main/splash-window-manager.ts:34`, `src/main/splash-window-manager.ts:35`
- Electron セッションで CSP と権限拒否のハンドラが登録されています: `src/main/session-hardening.ts:35`, `src/main/session-hardening.ts:83`, `src/main/session-hardening.ts:98`
- 外部 URL は IPC 経由で `http:` / `https:` のみに制限されています: `src/main/ipc-security.ts:31`, `src/main/ipc-security.ts:33`, `src/main/ipc-security.ts:69`
- Renderer から backend への `X-Auth-Token` は main process の network layer で注入され、renderer に直接露出していません: `src/main/backend-auth-injection.ts:18`, `src/main/backend-auth-injection.ts:61`, `src/main/backend-auth-injection.ts:78`
- FastAPI は docs/openapi をデフォルト無効化し、認証、CORS、リクエストサイズ、Host 検証を組み込んでいます: `src/python/app.py:171`, `src/python/app.py:240`, `src/python/app.py:390`, `src/python/app.py:430`
- Gemini API key は HTTP レスポンスでマスクされ、POSIX では設定ファイル権限を `0600` に寄せています: `src/python/api/settings.py:19`, `src/python/api/settings.py:35`, `src/python/quantum_calc/settings_manager.py:15`, `src/python/quantum_calc/settings_manager.py:247`
- React Markdown は `rehype-raw` を使っておらず、公式仕様上の既定 URL transform により `javascript:` URL はブロックされます: `src/web/components/ChatMessage.tsx:74`, `src/web/components/ChatMessage.tsx:76`

## High Severity

### SEC-001: Electron 37 系が公式サポート終了済み

Rule ID: ELECTRON-SUPPLY-001

Severity: High

Location:

- `package.json:68`
- `package-lock.json:5962`
- `package-lock.json:5963`

Evidence:

```json
"electron": "^37.2.4"
```

```json
"node_modules/electron": {
  "version": "37.7.0"
}
```

Impact:

Electron は Chromium と Node.js を同梱するため、EOL の Electron 37 系を使い続けると、Chromium/Node.js/Electron 本体のセキュリティ修正を受けられず、XSS や renderer compromise の影響が大きくなります。

Details:

- 公式 Electron release schedule では Electron 37.0.0 の EOL は 2026-01-13 です。
- 2026-06-09 時点の公式 stable releases では Electron 42.3.3 が最新として表示され、公式サポートポリシーは最新 3 stable major versions のみを対象にしています。

Fix:

- Electron をサポート中の major に上げる。現時点では 42 系を第一候補にし、互換性リスクが大きい場合でも 40/41/42 のサポート対象内に移行する。
- 併せて `electron-builder`, `electron-vite`, native module rebuild、packaging、`npm run test:main`, `npm run test:web`, `npm run build:electron` を検証する。

Mitigation:

- すぐに major upgrade できない場合でも、Electron 37 内の最新 patch に留まっていることを確認しつつ、移行 issue/期限を明示する。ただし EOL 後の security fix は期待できません。

False positive notes:

- `package.json` は `^37.2.4` ですが、lockfile は `37.7.0` に固定しています。いずれも 37 系であり、公式 EOL には該当します。

## Medium Severity

### SEC-002: CI がロック済み Python 環境を直接 `pip install` で変更している

Rule ID: SUPPLY-LOCK-001

Severity: Medium

Location:

- `.github/workflows/ci.yml:51`
- `.github/workflows/ci.yml:55`
- `.github/workflows/ci.yml:108`
- `.github/environment.yml:29`

Evidence:

```yaml
conda-lock install --name pyscf-env .github/pyscf-env.conda-lock.yml
```

```yaml
pip install "setuptools<75"
```

```yaml
- setuptools<75
```

Impact:

CI の build step が lockfile で作った環境を後から unpinned range の `pip install` で変更するため、再現性と供給網監査性が落ちます。PyPI 側の解決結果が変わると、同じコミットでも異なる依存物で検証・成果物作成が行われる可能性があります。

Fix:

- `.github/workflows/ci.yml:108` の直接 `pip install` を削除する。
- 必要な setuptools 制約はすでに `.github/environment.yml:29` にあるため、変更が必要なら `.github/environment.yml` を更新し、`npm run conda-lock:generate` で lockfile を再生成する。
- どうしても CI 内で追加導入が必要な場合は、完全固定バージョンとハッシュ検証、または conda-lock による管理に寄せる。

Mitigation:

- CI の build artifact を配布対象にしない運用なら影響は限定的です。ただし同じ workflow 内で成果物 upload があるため、ロック外依存の混入は避ける方が安全です。

False positive notes:

- Release workflow はこの直接 `pip install` を含まず、`npm run build` / `npm run build:linux` 経由です。ただし CI workflow では lockfile 後に明示的な環境変更が存在します。

### SEC-003: 依存関係の脆弱性スキャンが CI に見えない

Rule ID: SUPPLY-AUDIT-001

Severity: Medium

Location:

- `.github/workflows/ci.yml:57`
- `.github/workflows/ci.yml:89`
- `.github/workflows/ci.yml:170`
- `.github/workflows/ci.yml:179`

Evidence:

```yaml
npm ci --legacy-peer-deps
```

```yaml
python -m pytest tests/ -v
```

```yaml
npx tsc --noEmit
```

```yaml
npm run format:check
```

Impact:

Lockfile はありますが、既知脆弱性の検出ステップが CI に見えません。Electron の EOL 検出、npm パッケージ、FastAPI/Starlette/Uvicorn などの advisory 追跡が手作業頼みになり、古い依存が長く残る可能性があります。

Fix:

- Node: `npm audit --audit-level=high`、GitHub dependency-review、または OSV-Scanner を CI に追加する。
- Python/conda: `pip-audit` は conda 管理パッケージの扱いに注意が必要なため、OSV-Scanner や GitHub dependency graph/Dependabot の併用を検討する。
- Electron は EOL チェックを追加する。例: lockfile の `electron` major と公式 schedule を定期的に比較する運用、または Dependabot/Renovate で major update を issue 化する。

Mitigation:

- ネットワークが不安定な CI では、audit job を build/test と分離し、失敗ポリシーを段階的に導入する。

False positive notes:

- この監査環境では `npm audit --audit-level=high --omit=dev --json` が `registry.npmjs.org` 名前解決不可で失敗しました。したがって、現時点の advisory 有無は未確認です。

## Low Severity

### SEC-004: write request schema が未知フィールドを明示拒否していない

Rule ID: FASTAPI-VALID-001

Severity: Low

Location:

- `src/api-spec/openapi.yaml:1776`
- `src/api-spec/openapi.yaml:1798`
- `src/api-spec/openapi.yaml:3161`
- `src/api-spec/openapi.yaml:3897`
- `src/api-spec/openapi.yaml:4031`
- `src/api-spec/openapi.yaml:4041`
- `src/python/generated_models.py:64`
- `src/python/generated_models.py:86`
- `src/python/generated_models.py:1629`
- `src/python/generated_models.py:1703`
- `src/python/generated_models.py:1709`

Evidence:

```yaml
QuantumCalculationRequest:
  oneOf:
```

```yaml
SettingsUpdateRequest:
  $ref: '#/components/schemas/AppSettings'
```

```python
class AgentChatRequest(BaseModel):
```

```python
class CreateChatSessionRequest(BaseModel):
```

Impact:

現在は Pydantic のモデル化後に `model_dump` を使っており、未知フィールドが直接 DB 更新に流れる mass assignment は確認していません。ただし API 契約が未知フィールドを明示拒否しないため、クライアントの誤送信や将来の実装変更が見落とされやすくなります。

Fix:

- OpenAPI の write schema に `additionalProperties: false` を追加し、`npm run codegen` で Python/TypeScript 生成物を更新する。
- 生成モデルで `extra='forbid'` 相当になることを確認する。
- まず `SettingsUpdateRequest`, `AgentChatRequest`, `CreateChatSessionRequest`, `UpdateChatSessionRequest`, `CalculationUpdateRequest` など、状態変更系から適用する。

Mitigation:

- 現状維持の場合でも、各 handler で Pydantic モデルから dump した値だけを使う方針を維持し、raw dict を直接永続化しない。

False positive notes:

- 現行実装では未知フィールドは直接利用されていないため、これは即時 exploit ではなく hardening finding です。

## Notable Non-Findings

- `dangerouslySetInnerHTML`, `innerHTML`, `eval`, `new Function`, `postMessage` の危険な実使用は確認できませんでした。
- `localStorage` は `pyscf_setup_completed` フラグのみで、認証 token / session / API key の保存は確認できませんでした: `src/web/App.tsx:79`, `src/web/App.tsx:113`
- SQL は sqlite の parameterized query を使っています: `src/python/database/chat_history.py:342`, `src/python/database/chat_history.py:414`, `src/python/database/chat_history.py:450`
- 計算 ID の path traversal は `CalculationRepository.resolve_calculation_path` 経由で防御され、回帰テストもあります: `src/python/services/calculation_service_context.py:73`, `src/python/tests/integration/test_api_endpoints/test_quantum_api.py:552`
- Swagger UI は development-only の設計で、SRI と CSP が入っています: `src/python/api/swagger_ui.py:4`, `src/python/api/swagger_ui.py:30`, `src/python/api/swagger_ui.py:48`, `src/python/api/swagger_ui.py:91`, `src/python/api/swagger_ui.py:99`

## Verification Notes

実行済み:

- `rg` による高リスクパターン検索
- `package-lock.json` から Electron 実解決バージョン確認
- `npm audit --audit-level=high --omit=dev --json`

未完了:

- `npm audit` はネットワーク制限により失敗しました: `getaddrinfo ENOTFOUND registry.npmjs.org`
- Python/conda の advisory scan は実行していません。conda lockfile 全体を扱える scanner の選定が必要です。
- 実修正は行っていません。このレポートは監査結果のみです。

## Sources

- Electron Security checklist: https://www.electronjs.org/docs/latest/tutorial/security
- Electron release/support policy: https://www.electronjs.org/docs/latest/tutorial/electron-timelines
- Electron release schedule: https://releases.electronjs.org/schedule
- Electron stable releases: https://releases.electronjs.org/releases/stable
- react-markdown README/security behavior: https://github.com/remarkjs/react-markdown/blob/main/readme.md
