# 設計

## 概要
PySCF ドキュメントに記載されている主要な基底関数・XC 汎函数を `supported_parameters` に追加し、UI がそのまま拡張された選択肢を表示できるようにする。

## データフロー
- `src/python/quantum_calc/supported_parameters.py` にあるカテゴリ別リストを拡張。
- `get_all_supported_parameters()` が API に返却。
- 計算条件 UI がカテゴリごとに選択肢を描画。

## 変更点
- `get_supported_basis_functions()` に追加候補の基底関数を追記。
- `get_supported_exchange_correlation()` に追加候補の XC 汎函数を追記。
- UI 側の構造変更は不要（カテゴリ名は既存/追加分をそのまま利用）。

## 追加する項目
- 要件の「追加候補（提案）」に記載したリストをそのまま反映する。
- CAM-B3LYP の表記で統一する。
- IR スペクトルのスケールファクターは追加せず、警告付きフォールバックを許容する。
- PySCF/libxc の対応状況に依存するため、動的検証の追加は行わない。

## テスト方針
- 手動: 計算条件画面で追加候補がカテゴリ別に表示されること。
- 既存テストにリスト固定の期待値があれば更新。

### コマンド
- `npm run format:check`
- `npm run test:build`
- `npm run test:python-build`
- `/Users/goodapple/miniforge3/envs/pyscf-env/bin/python -m pytest src/python/tests -v`
