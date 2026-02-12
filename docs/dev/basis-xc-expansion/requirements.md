# 要件

## 目的
PySCF のドキュメントに基づき、アプリで選択できる基底関数と交換相関（XC）汎函数を主要なものへ拡張する。

## ユーザーストーリー
- ユーザーとして、用途に合わせてより多くの基底関数・XC 汎函数を選べる。

## 対象外
- 計算アルゴリズムや PySCF 本体の変更。
- 既存のデフォルト設定の変更。
- UI の再設計（選択肢の拡張のみ）。

## 追加候補（提案）
### 基底関数
- Minimal
  - STO-3G
  - STO-6G
  - MINAO
- Pople（分割原子価）
  - 3-21G
  - 4-31G
  - 6-31G
  - 6-31G(d)
  - 6-31G(d,p)
  - 6-31+G(d)
  - 6-31+G(d,p)
  - 6-31++G(d,p)
  - 6-311G
  - 6-311G(d)
  - 6-311G(d,p)
  - 6-311+G(d)
  - 6-311+G(d,p)
  - 6-311++G(d,p)
- Dunning（相関一貫）
  - cc-pVDZ
  - cc-pVTZ
  - cc-pVQZ
  - cc-pV5Z
  - aug-cc-pVDZ
  - aug-cc-pVTZ
  - aug-cc-pVQZ
  - aug-cc-pV5Z
- def2
  - def2-SVP
  - def2-SVPD
  - def2-TZVP
  - def2-TZVPP
  - def2-TZVPD
  - def2-QZVP
  - def2-QZVPP
  - def2-QZVPD

### 交換相関（XC）汎函数
- LDA
  - SVWN
- GGA
  - PBE
  - BLYP
  - BP86
  - PW91
  - revPBE
  - PBEsol
  - RPBE
  - OLYP
  - BPBE
- Meta-GGA
  - M06-L
  - TPSS
  - revTPSS
  - SCAN
  - rSCAN
  - r2SCAN
- Hybrid (GGA)
  - B3LYP
  - B3PW91
  - PBE0
  - B3P86
  - O3LYP
- Hybrid (Meta-GGA)
  - M06
  - M06-2X
  - TPSSh
- Range-separated Hybrid
  - CAM-B3LYP
  - wB97XD

※ 既存の選択肢は維持する。

## 受入条件
- supported-parameters API が拡張後の基底関数・XC 汎函数を返す。
- 計算条件 UI に追加候補がカテゴリ別に表示される。
- 既存の選択肢は削除されない。
 - 実環境の PySCF/libxc の対応状況により、選択した項目が実行時にエラーとなる可能性は許容する（動的検証は今回行わない）。

## 確認事項
- 上記の追加候補で進めてよい（承認済み）。
- CAM-B3LYP の表記で統一する（承認済み）。
- IR スペクトルのスケールファクターは警告付きフォールバックを許容（承認済み）。
