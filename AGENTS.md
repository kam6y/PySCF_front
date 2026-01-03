あなたのパートナーは日本人です。なので報告は日本語でお願いします。開発は以下の手順で進めてください。

1. docs/dev/{feature_name}/{requirements,design}.md を作成する。質問があれば私に聞いてください。最終的に私の承認をとってください。
2. 開発を実施。
3. lint, build, test を pass する。
   - lint: npm run format:check（失敗したら npm run format → npm run format:check）
   - build: npm run test:build
   - test: npm run test:python-build と pytest を必ず実行
     - pytest は conda 環境の Python を固定で使用する
       - /Users/goodapple/miniforge3/envs/pyscf-env/bin/python -m pytest src/python/tests -v
4. Bash で "my-review" コマンドの標準入力にメッセージを送信するとレビューしてもらえます。
   - 推奨: npm run my-review（デフォルト 900 秒タイムアウト）
   - my-review は codex review を使い、未コミット差分（--uncommitted）を自動レビューする
   - タイムアウトを延長したい場合:
     - MY_REVIEW_TIMEOUT=1800 npm run my-review
     - npm run my-review -- --timeout 1800
   - 実行例（パイプ）:
     - cat <<'EOF' | npm run my-review
       docs/dev/{feature_name}/requirements.md と design.md を更新し、<変更内容の要約> を実装しました。
       主要差分: <主要差分の要約>
       重点確認: 設計妥当性、潜在バグ、足りないテスト
       EOF
   docs/dev/{feature_name}/{requirements,design}.md と開発内容を伝え、main との差分の品質・設計妥当性・潜在バグ・不足テストなどを指摘してもらってください。
   レビューは時間がかかるのでデフォルトは 900 秒です。必要なら上記の方法で延長してください。
5. 妥当な改善点であれば実施し、もう一度 4. を実施する。"review→改善"の繰り返しは最大2回まで。
