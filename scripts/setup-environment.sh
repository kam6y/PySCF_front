#!/bin/bash

# PySCF Native App - 自動環境構築スクリプト
# このスクリプトは、開発に必要なconda環境を自動的にセットアップします

set -e  # エラー時に終了

# カラー出力用の定数
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ログ出力関数
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# プロジェクトルートディレクトリの取得
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

log_info "PySCF Native App 環境構築を開始します..."
log_info "プロジェクトルート: $PROJECT_ROOT"

# conda の存在確認
if ! command -v conda &> /dev/null; then
    log_error "conda が見つかりません。"
    log_info "conda のインストール方法:"
    log_info "1. Miniforge (推奨): https://github.com/conda-forge/miniforge"
    log_info "2. Miniconda: https://docs.conda.io/en/latest/miniconda.html"
    log_info "3. Anaconda: https://www.anaconda.com/products/distribution"
    exit 1
fi

log_success "conda が見つかりました: $(conda --version)"

# conda の初期化確認
CONDA_BASE=$(conda info --base 2>/dev/null)
if [[ -n "$CONDA_BASE" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    log_info "conda環境を読み込みました: $CONDA_BASE"
else
    log_warning "conda初期化スクリプトが見つかりません"
    log_info "手動でcondaを初期化してください: conda init"
fi

# 既存環境の確認
ENV_NAME="pyscf-env"
if conda env list | grep -q "^$ENV_NAME "; then
    log_warning "環境 '$ENV_NAME' は既に存在します"
    read -p "既存の環境を削除して再作成しますか？ (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        log_info "既存環境を削除しています..."
        conda env remove -n "$ENV_NAME" -y
        log_success "既存環境を削除しました"
    else
        log_info "既存環境をそのまま使用します"
        conda activate "$ENV_NAME"
        log_success "環境 '$ENV_NAME' をアクティベートしました"
        exit 0
    fi
fi

# environment.yml ファイルの確認
ENV_FILE="$PROJECT_ROOT/.github/environment.yml"
if [[ ! -f "$ENV_FILE" ]]; then
    log_error "environment.yml が見つかりません: $ENV_FILE"
    exit 1
fi

log_info "environment.yml を使用して環境を作成します: $ENV_FILE"

# conda環境の作成
log_info "conda環境 '$ENV_NAME' を作成中..."
if conda env create -f "$ENV_FILE"; then
    log_success "conda環境の作成が完了しました"
else
    log_error "conda環境の作成に失敗しました"
    exit 1
fi

# 環境のアクティベート
log_info "環境をアクティベートしています..."
conda activate "$ENV_NAME"

# 環境の検証
log_info "環境の検証を実行中..."
# conda環境のPythonを直接使用して、新しく作成した環境を検証
"$(conda info --base)/envs/$ENV_NAME/bin/python" "$PROJECT_ROOT/scripts/verify-environment.py"

if [[ $? -eq 0 ]]; then
    log_success "環境の検証が成功しました"
else
    log_error "環境の検証に失敗しました"
    exit 1
fi

# Node.js依存関係のインストール確認
log_info "Node.js依存関係を確認中..."
cd "$PROJECT_ROOT"
if [[ ! -d "node_modules" ]]; then
    log_info "Node.js依存関係をインストールしています..."
    npm install
    log_success "Node.js依存関係のインストールが完了しました"
else
    log_info "Node.js依存関係は既にインストールされています"
fi

# 成功メッセージ
log_success "環境構築が完了しました！"
echo
log_info "次のステップ:"
log_info "1. 環境をアクティベート: conda activate $ENV_NAME"
log_info "2. 開発サーバーを起動: npm run dev"
echo
log_info "環境変数でconda環境パスを指定したい場合:"
log_info "export CONDA_ENV_PATH=\$(conda info --base)/envs/$ENV_NAME"