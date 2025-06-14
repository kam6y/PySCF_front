# Database Docker Configuration

このディレクトリは、MySQLのカスタムDockerfileが必要になった場合に使用します。

## 現在の状況

現在は公式MySQL 8.0イメージを直接使用しているため、カスタムDockerfileは不要です。

## 初期化について

データベースの初期化は以下で自動実行されます：
- `../database/schema/001_initial_schema.sql` → `/docker-entrypoint-initdb.d/`

## カスタムDockerfileが必要になる場合

以下のような要件がある場合、このディレクトリにDockerfileを作成します：

- カスタムMySQL設定ファイル
- 追加的なMySQLプラグインのインストール
- 特別なユーザー権限設定
- カスタム初期化スクリプト

## 使用例

```dockerfile
# 将来的にカスタマイズが必要な場合の例
FROM mysql:8.0

# カスタム設定ファイル
COPY custom-my.cnf /etc/mysql/conf.d/

# 追加的な初期化スクリプト
COPY init-scripts/ /docker-entrypoint-initdb.d/
```