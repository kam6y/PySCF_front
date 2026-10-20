# macOS 配布向けの署名と公証 (Code Signing / Notarization)

macOS 向け DMG を Gatekeeper の警告なしで配布するには、

1. **Developer ID Application** 証明書で署名 (Hardened Runtime 有効)
2. Apple の **公証 (Notarization)** を通し、チケットをアプリにステープル

の両方が必要です。このリポジトリでは `electron-builder` がこの 2 つを自動で行い、
`scripts/verify-mac-signing.sh` が結果を検証します。設定は以下にあります。

| 役割 | 場所 |
| --- | --- |
| 署名・公証の設定 | `package.json` → `build.mac` (`hardenedRuntime`, `entitlements*`, `notarize`) |
| Hardened Runtime の例外 (本体) | `src/assets/mac/entitlements.mac.plist` |
| Hardened Runtime の例外 (同梱 Python など) | `src/assets/mac/entitlements.mac.inherit.plist` |
| 検証スクリプト | `scripts/verify-mac-signing.sh` (`npm run verify:mac-signing`) |
| CI (タグ push で実行) | `.github/workflows/release.yml` |

> **注意**: `Apple Distribution` 証明書は App Store 提出専用で、DMG 配布や公証には使えません。
> `Apple Development` 証明書は手元での動作確認には使えますが、配布物には使えません。
> 配布には必ず **Developer ID Application** 証明書が必要です。

## 1. 一度だけ行う準備 (Apple 側)

### 1-1. Developer ID Application 証明書を作成してインストールする

Developer ID 証明書は Apple Developer Program の **Account Holder** しか作成できません。

**Xcode から作る場合** (推奨):

1. Xcode → Settings… → Accounts → 対象の Team を選択 → *Manage Certificates…*
2. 左下の **+** → **Developer ID Application**
3. 作成されると自動でログインキーチェーンに入ります。

**Web から作る場合**:

1. キーチェーンアクセス → 証明書アシスタント → *認証局に証明書を要求…* で CSR (`.certSigningRequest`) を作成
2. <https://developer.apple.com/account/resources/certificates/add> → **Developer ID Application** → CSR をアップロード
3. ダウンロードした `.cer` をダブルクリックしてログインキーチェーンに追加

確認:

```bash
security find-identity -v -p codesigning
#  ... "Developer ID Application: <名前> (<TEAM_ID>)" が表示されれば OK
```

### 1-2. 公証用の認証情報を用意する

どちらか一方を用意します。

**A. Apple ID + App 用パスワード (簡単)**

1. <https://account.apple.com/account/manage> → *サインインとセキュリティ* → *App 用パスワード* で生成
2. Team ID を確認: <https://developer.apple.com/account> → *Membership details*
   (`security find-identity` に表示される括弧内の英数字と同じです)

**B. App Store Connect API キー (electron-builder 推奨)**

1. <https://appstoreconnect.apple.com/access/integrations/api> → *Team Keys* → **+**
   (ロールは *Developer* 以上)
2. `AuthKey_XXXXXXXXXX.p8` をダウンロード (再ダウンロード不可)。**Key ID** と **Issuer ID** も控える

## 2. 手元で署名付きビルドを試す

```bash
# 署名のみ (Developer ID 証明書がキーチェーンにあれば自動で選択されます)
npm run package

# 公証まで行う場合は環境変数を付けて electron-builder を実行
export APPLE_ID="you@example.com"
export APPLE_APP_SPECIFIC_PASSWORD="xxxx-xxxx-xxxx-xxxx"
export APPLE_TEAM_ID="XXXXXXXXXX"
npm run build && npx electron-builder --mac && npm run verify:asar-artifacts

# 結果を検証 (署名・Hardened Runtime・公証チケット・Gatekeeper)
npm run verify:mac-signing
```

- 認証情報を毎回 export したくない場合は notarytool のキーチェーンプロファイルを使えます:
  `xcrun notarytool store-credentials pyscf-notary --apple-id ... --team-id ... --password ...`
  を一度実行し、以後は `APPLE_KEYCHAIN_PROFILE=pyscf-notary` を付けて electron-builder を実行します。
- Developer ID 証明書がまだ無い環境では、電子署名の流れだけを
  `npm run verify:mac-signing -- --allow-development` で確認できます (配布不可)。
- 公証は 500 MB 超の DMG をアップロードするため、数分〜十数分かかります。
  失敗時は `xcrun notarytool log <submission-id> --keychain-profile pyscf-notary` で理由を確認できます。

## 3. GitHub Actions の Secrets を登録する

`v*` タグを push すると `.github/workflows/release.yml` が macOS ランナー上で署名・公証・検証まで行います。
以下の Secrets を **Settings → Secrets and variables → Actions** に登録してください
(`gh secret set <NAME>` でも登録できます)。

### 3-1. 署名用 (必須)

| Secret | 内容 |
| --- | --- |
| `CSC_LINK` | Developer ID Application 証明書 (+秘密鍵) を書き出した `.p12` を base64 にしたもの |
| `CSC_KEY_PASSWORD` | `.p12` 書き出し時に設定したパスワード |

`.p12` の作り方: キーチェーンアクセス → *自分の証明書* → **Developer ID Application: …** を右クリック → *書き出す…* → フォーマット `.p12`。

```bash
base64 -i DeveloperID.p12 | gh secret set CSC_LINK
gh secret set CSC_KEY_PASSWORD      # 対話的に入力
```

### 3-2. 公証用 (A か B のどちらか)

| 方式 | Secret |
| --- | --- |
| A. Apple ID | `APPLE_ID`, `APPLE_APP_SPECIFIC_PASSWORD`, `APPLE_TEAM_ID` |
| B. API キー | `APPLE_API_KEY_CONTENT` (`.p8` ファイルの中身そのまま), `APPLE_API_KEY_ID`, `APPLE_API_ISSUER` |

```bash
# A
gh secret set APPLE_ID
gh secret set APPLE_APP_SPECIFIC_PASSWORD
gh secret set APPLE_TEAM_ID
# B
gh secret set APPLE_API_KEY_CONTENT < AuthKey_XXXXXXXXXX.p8
gh secret set APPLE_API_KEY_ID
gh secret set APPLE_API_ISSUER
```

Secrets が足りない場合、ワークフローは「Check macOS code-signing secrets」ステップで
不足している名前を表示して失敗します (署名なしの DMG が公開されるのを防ぐため)。

## 4. リリース手順

```bash
# main に変更を取り込んだ上で (従来どおり)
npm run release -- 0.1.8-beta
```

タグ `v0.1.8-beta` が push されると Release ワークフローが

1. Linux (AppImage) と macOS (DMG) をビルド
2. macOS は Developer ID で署名 → 公証 → ステープル → `scripts/verify-mac-signing.sh` で検証
3. GitHub Release に成果物を添付

します。利用者は DMG を開き、アプリを Applications にドラッグするだけで起動できます
(右クリック →「開く」や `xattr -d com.apple.quarantine` は不要になります)。

## 5. 仕組みのメモ

- `electron-builder` は `.app` 内の全 Mach-O (Electron 本体・ヘルパー・`Resources/conda_env` 以下の
  `python3.12` や `.so`/`.dylib` 約 900 個) を個別に署名します。`codesign` は 1 ファイルごとに Apple の
  タイムスタンプサーバーへ問い合わせるため、1 ファイルあたり 0.2〜0.4 秒かかります。
- 署名ツールは「テキストでないファイル」をすべて署名対象にするため、何もしないと terminfo や zoneinfo、
  フォント、画像など Mach-O ではない約 4,000 ファイルにも署名が走り、30 分以上かかります。
  `package.json` の `build.mac.signIgnore` で「Mach-O になり得ないパス・拡張子」だけを除外し、
  署名を数分に抑えています。Mach-O を誤って除外していないことは `scripts/verify-mac-signing.sh` が検出します。
- conda 環境内のディレクトリへのシンボリックリンク (`lib/python3.1 -> python3.12`、`lib/terminfo`、
  `lib/icu/current`) も署名ツールが辿るため、放置すると site-packages 全体が二重に署名されます。
  これらの別名パスも `signIgnore` で除外しています (実体のパスは通常どおり署名されます)。
- `codesign` が `The timestamp service is not available.` で失敗することがあります。これは Apple の
  タイムスタンプサーバー側の一時的な障害で、electron-builder は 1 ファイルあたり 3 回まで再試行します。
  それでも失敗した場合は設定の問題ではないので、しばらく待ってからワークフロー (またはローカルのビルド)
  を再実行してください。
- 実行時に不要な `*.o` (Mach-O オブジェクトファイル) は `extraResources` のフィルタで同梱対象から外しています。
  一方 `etc/conda/test-files/` などは不要に見えても削除してはいけません。初回起動時の `conda-unpack` は
  自身の置換リストにあるファイルが 1 つでも欠けると `FileNotFoundError` で失敗し、アプリが起動できなくなります
  (同梱内容を減らしたい場合は conda-pack の前に環境側で削除してください)。
- 初回起動時にアプリは `conda_env` を `~/Library/Application Support/Pyscf_front/conda_env` へコピーし
  `conda-unpack` を実行します。`conda-unpack` はプレフィックスを書き換えたバイナリを **ad-hoc で再署名**するため、
  同梱 Python がそれらを読み込めるように `entitlements.mac.inherit.plist` で
  `com.apple.security.cs.disable-library-validation` を付けています。
- electron-builder は **26.16.1 以上** が必要です。それより前のバージョンは CI 用の一時キーチェーンを作るときに
  `.p12` のパスワードをキーチェーンのパスワードとして `security set-key-partition-list` に渡すバグがあり
  (electron-builder #10066)、macOS 26 以降のランナーでは
  `The user name or passphrase you entered is not correct.` で署名前に失敗します。
- DMG 自体は署名・公証していません (electron-builder の既定)。Gatekeeper は DMG 内の `.app` に
  ステープルされたチケットを見るため、これで問題なく起動できます。
