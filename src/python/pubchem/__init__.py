"""PubChem APIから分子構造を取得するためのクライアントモジュールです。"""

# メインのクライアントとデータクラスをパッケージから直接公開します
from .client import PubChemClient, PubChemError, CompoundData

# parserモジュールをインポートしてアクセスできるようにします
from . import parser

__version__ = "0.2.0"

# __all__ を更新して、新しい正しいエクスポートを反映させます
__all__ = [
    "PubChemClient",
    "PubChemError",
    "CompoundData",
    "parser"
]