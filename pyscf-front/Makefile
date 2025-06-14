.PHONY: help setup-backend setup-frontend setup test run-backend run-frontend docker-up docker-down clean

help:
	@echo "使用可能なコマンド:"
	@echo "  make setup          - 開発環境の完全セットアップ"
	@echo "  make setup-backend  - バックエンド環境のセットアップ"
	@echo "  make setup-frontend - フロントエンド環境のセットアップ"
	@echo "  make test          - 全テストの実行"
	@echo "  make run-backend   - バックエンドサーバーの起動"
	@echo "  make run-frontend  - フロントエンドの起動"
	@echo "  make docker-up     - Docker環境の起動"
	@echo "  make docker-down   - Docker環境の停止"
	@echo "  make clean         - キャッシュファイルの削除"

setup: setup-backend setup-frontend
	@echo "開発環境のセットアップが完了しました"

setup-backend:
	cd backend && python -m venv venv
	cd backend && ./venv/bin/pip install --upgrade pip
	cd backend && ./venv/bin/pip install -r requirements.txt
	cd backend && ./venv/bin/pip install -e .

setup-frontend:
	cd frontend && flutter pub get
	cd frontend && dart run build_runner build

test:
	cd backend && ./venv/bin/pytest
	cd frontend && flutter test

run-backend:
	cd backend && ./venv/bin/python -m uvicorn src.main:app --reload

run-frontend:
	cd frontend && flutter run -d macos  # または windows, linux

docker-up:
	docker-compose -f docker-compose.yml up -d

docker-up-dev:
	docker-compose -f docker-compose.yml -f docker-compose.dev.yml up -d

docker-down:
	docker-compose -f docker-compose.yml down

docker-build:
	docker-compose -f docker-compose.yml build

clean:
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete
	find . -type d -name ".pytest_cache" -delete
	cd frontend && flutter clean