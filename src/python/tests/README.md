# PySCF Front Backend Test Suite

包括的なテスト戦略に基づいた、PySCF Frontバックエンドのテストスイートです。

## 📁 ディレクトリ構造

```
tests/
├── README.md                    # このファイル
├── conftest.py                  # 共通フィクスチャ定義
├── test_fixtures.py            # テスト基盤の動作検証
├── data/                       # テストデータとモック
│   ├── mock_pubchem_response.json
│   ├── sample_h2.xyz
│   └── sample_water.xyz
├── unit/                       # 単体テスト（隔離環境）
│   ├── test_services/         # サービス層のテスト
│   ├── test_quantum_calc/     # 量子化学計算モジュールのテスト
│   └── test_agent/            # AIエージェントのテスト
└── integration/                # 結合テスト（コンポーネント連携）
    └── test_api_endpoints/    # APIエンドポイントのテスト
```

## 🚀 テストの実行

### 基本的な実行

```bash
# conda環境をアクティベート
conda activate pyscf-env

# src/pythonディレクトリに移動
cd src/python

# すべてのテストを実行
pytest tests/ -v

# 特定のテストファイルを実行
pytest tests/test_fixtures.py -v

# 特定のテスト関数を実行
pytest tests/test_fixtures.py::test_app_fixture -v

# マーカーでフィルタリング
pytest tests/ -m unit -v           # 単体テストのみ
pytest tests/ -m integration -v    # 結合テストのみ
pytest tests/ -m "not slow" -v     # 低速テストを除外
```

### カバレッジレポート付き実行

```bash
# ターミナルにカバレッジサマリーを表示
pytest tests/ --cov=. --cov-report=term-missing

# HTMLカバレッジレポートを生成（htmlcov/index.htmlをブラウザで開く）
pytest tests/ --cov=. --cov-report=html

# 両方を同時に生成
pytest tests/ --cov=. --cov-report=term-missing --cov-report=html
```

### 便利なオプション

```bash
# 最初の失敗で停止
pytest tests/ -x

# 失敗したテストのみ再実行
pytest tests/ --lf

# より詳細な出力
pytest tests/ -vv

# テストの実行時間を表示
pytest tests/ --durations=10

# 並列実行（pytest-xdistが必要）
pytest tests/ -n auto
```

## 🧪 利用可能なフィクスチャ

`conftest.py`で定義された主要なフィクスチャ：

### アプリケーション関連

- **app**: テスト用Flaskアプリケーションインスタンス（session scope）
- **client**: HTTPリクエスト用テストクライアント（function scope）
- **socketio_client**: WebSocket通信用テストクライアント（function scope）
- **runner**: Flask CLIコマンドテスト用ランナー

### テストデータ

- **sample_h2_xyz**: 水素分子のXYZ座標
- **sample_water_xyz**: 水分子のXYZ座標
- **valid_dft_params**: 有効なDFT計算パラメータ
- **valid_hf_params**: 有効なHartree-Fock計算パラメータ

### ヘルパー

- **dummy_executor**: 非同期処理の同期テスト用Executorクラスのインスタンス

### 使用例

```python
def test_api_endpoint(client):
    """HTTPリクエストのテスト例"""
    response = client.get('/api/quantum/calculations')
    assert response.status_code == 200
    data = response.get_json()
    assert data['success'] is True


def test_websocket_event(socketio_client):
    """WebSocketのテスト例"""
    socketio_client.emit('join_calculation', {'calculation_id': 'test-123'})
    received = socketio_client.get_received()
    assert len(received) > 0


def test_with_sample_data(sample_h2_xyz, valid_dft_params):
    """サンプルデータを使用したテスト例"""
    assert 'H' in sample_h2_xyz
    assert valid_dft_params['calculation_method'] == 'DFT'
```

## 📝 テストの書き方ガイドライン

### テストの命名規則

- ファイル名: `test_*.py`
- テスト関数: `test_*`
- テストクラス: `Test*`

### テスト構造（Given-When-Then）

```python
def test_calculation_success(client, valid_dft_params):
    """
    GIVEN 有効な計算パラメータ
    WHEN 計算開始APIが呼び出される
    THEN 202 Acceptedが返され、計算インスタンスが作成される
    """
    # ARRANGE (Given)
    payload = valid_dft_params

    # ACT (When)
    response = client.post('/api/quantum/calculate', json=payload)

    # ASSERT (Then)
    assert response.status_code == 202
    data = response.get_json()
    assert data['success'] is True
    assert 'calculation' in data['data']
```

### マーカーの使用

```python
import pytest

@pytest.mark.unit
def test_isolated_logic():
    """単体テストのマーカー例"""
    pass

@pytest.mark.integration
def test_component_interaction():
    """結合テストのマーカー例"""
    pass

@pytest.mark.slow
def test_heavy_computation():
    """低速テストのマーカー例"""
    pass

@pytest.mark.parametrize("input,expected", [
    (1, 2),
    (2, 4),
    (3, 6)
])
def test_with_parameters(input, expected):
    """パラメータ化テストの例"""
    assert input * 2 == expected
```

## 🎯 テスト戦略

### 単体テスト（Unit Tests）

- **目的**: 個々のコンポーネントのロジックを隔離して検証
- **特徴**:
  - 外部依存をモックで置き換え
  - 高速実行（< 1秒）
  - 高い再現性
- **対象**: サービス層、計算ロジック、ユーティリティ関数

### 結合テスト（Integration Tests）

- **目的**: 複数のコンポーネントの連携を検証
- **特徴**:
  - 実際のFlaskアプリケーションを使用
  - APIエンドポイントのE2Eテスト
  - WebSocket通信のテスト
- **対象**: APIエンドポイント、WebSocketハンドラ、ワークフロー全体

## 🔧 モッキング戦略

### pytest-mockの使用

```python
def test_with_mock(mocker):
    """外部APIをモック化する例"""
    # モックの設定
    mock_response = {"xyz": "H 0 0 0\nH 0 0 0.74"}
    mocker.patch('pubchem.client.PubChemClient.search', return_value=mock_response)

    # テスト実行
    from services.pubchem_service import PubChemService
    service = PubChemService()
    result = service.search_by_name("water")

    # 検証
    assert result == mock_response
```

### DummyExecutorによる同期テスト

```python
def test_async_workflow(mocker, client, valid_dft_params):
    """非同期処理を同期的にテストする例"""
    # ProcessPoolExecutorをDummyExecutorに置き換え
    from tests.conftest import DummyExecutor
    mocker.patch(
        'quantum_calc.process_manager.ProcessPoolExecutor',
        new=DummyExecutor
    )

    # 計算を開始（即座に完了）
    response = client.post('/api/quantum/calculate', json=valid_dft_params)
    assert response.status_code == 202

    calc_id = response.get_json()['data']['calculation']['id']

    # 結果を即座に取得できる
    result_response = client.get(f'/api/quantum/calculations/{calc_id}')
    assert result_response.status_code == 200
```

## 📊 カバレッジ目標

- **全体目標**: 80%以上
- **クリティカルパス**: 90%以上（計算ロジック、APIエンドポイント）
- **現在のカバレッジ**: HTMLレポート（`htmlcov/index.html`）で確認

## 🐛 トラブルシューティング

### ModuleNotFoundError

```bash
# pytest.iniがsrc/pythonディレクトリにあることを確認
cd src/python
pytest tests/
```

### WebSocket関連のエラー

```python
# SocketIOクライアントは自動的に接続を確立
def test_websocket(socketio_client):
    assert socketio_client.is_connected()  # 接続を確認
```

### カバレッジレポートが表示されない

```bash
# pytest-covがインストールされているか確認
conda list pytest-cov

# インストールされていない場合
conda install -c conda-forge pytest-cov
```

## 📚 参考資料

- [Pytest公式ドキュメント](https://docs.pytest.org/)
- [pytest-mock](https://pytest-mock.readthedocs.io/)
- [pytest-cov](https://pytest-cov.readthedocs.io/)
- [Flask Testing](https://flask.palletsprojects.com/en/latest/testing/)
- [Flask-SocketIO Testing](https://flask-socketio.readthedocs.io/en/latest/testing.html)

## 🔄 次のステップ

第I部（テスト基盤）が完了しました。次の実装項目：

1. **第II部**: 単体テストの実装
   - サービス層のテスト（PubChem, SMILES, 量子計算）
   - 計算モジュールのテスト（DFT, HF, MP2, CCSD等）
   - AIエージェントのテスト（Function Calling）

2. **第III部**: 結合テストの実装
   - APIエンドポイントのテスト
   - WebSocketハンドラのテスト
   - 計算ワークフロー全体のテスト

3. **第IV部**: CI/CDパイプラインへの統合
   - GitHub Actionsワークフローの更新
   - Codecov統合
   - 品質ゲートの設定
