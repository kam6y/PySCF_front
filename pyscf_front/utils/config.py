"""
設定管理モジュール
"""
import os
from pathlib import Path
from typing import Any, Dict, Optional
import yaml
from dotenv import load_dotenv

class Config:
    """設定管理クラス"""
    
    def __init__(self, config_file: Optional[str] = None):
        self.project_root = Path(__file__).parent.parent.parent
        self.config_file = config_file or self.project_root / "config.yaml"
        self.env_file = self.project_root / ".env"
        
        # 環境変数を読み込み
        if self.env_file.exists():
            load_dotenv(self.env_file)
        
        # 設定ファイルを読み込み
        self.config_data = self._load_config()
        
    def _load_config(self) -> Dict[str, Any]:
        """設定ファイルの読み込み"""
        default_config = {
            'database': {
                'url': 'sqlite:///pyscf_front.db',
                'echo': False
            },
            'application': {
                'debug': False,
                'log_level': 'INFO',
                'log_file': 'logs/pyscf_front.log'
            },
            'calculation': {
                'max_concurrent_jobs': 4,
                'default_method': 'B3LYP',
                'default_basis_set': '6-31G(d)'
            },
            'gui': {
                'theme': 'dark',
                'language': 'ja',
                'window_geometry': '1400x900'
            },
            'mcp_server': {
                'enabled': False,
                'port': 50053,
                'host': 'localhost'
            }
        }
        
        if isinstance(self.config_file, Path) and self.config_file.exists():
            try:
                with open(self.config_file, 'r', encoding='utf-8') as f:
                    file_config = yaml.safe_load(f)
                    if file_config:
                        default_config.update(file_config)
            except Exception as e:
                print(f"Warning: Could not load config file {self.config_file}: {e}")
        
        return default_config
    
    def get(self, key: str, default: Any = None) -> Any:
        """設定値の取得（環境変数を優先）"""
        # 環境変数をチェック（大文字でアンダースコア区切り）
        env_key = key.upper().replace('.', '_')
        env_value = os.getenv(env_key)
        if env_value is not None:
            return self._convert_env_value(env_value)
        
        # 設定ファイルから取得
        keys = key.split('.')
        value = self.config_data
        
        try:
            for k in keys:
                value = value[k]
            return value
        except (KeyError, TypeError):
            return default
    
    def _convert_env_value(self, value: str) -> Any:
        """環境変数の値を適切な型に変換"""
        if value.lower() in ('true', 'yes', '1', 'on'):
            return True
        elif value.lower() in ('false', 'no', '0', 'off'):
            return False
        elif value.isdigit():
            return int(value)
        else:
            try:
                return float(value)
            except ValueError:
                return value
    
    def set(self, key: str, value: Any):
        """設定値の更新"""
        keys = key.split('.')
        config = self.config_data
        
        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]
        
        config[keys[-1]] = value
    
    def save(self):
        """設定ファイルに保存"""
        try:
            if isinstance(self.config_file, Path):
                self.config_file.parent.mkdir(parents=True, exist_ok=True)
            with open(self.config_file, 'w', encoding='utf-8') as f:
                yaml.dump(self.config_data, f, default_flow_style=False, allow_unicode=True)
        except Exception as e:
            print(f"Error saving config: {e}")
    
    def get_database_url(self) -> str:
        """データベースURLの取得"""
        return self.get('database.url', 'sqlite:///pyscf_front.db')
    
    def is_debug(self) -> bool:
        """デバッグモードの確認"""
        return self.get('application.debug', False)
    
    def get_log_level(self) -> str:
        """ログレベルの取得"""
        return self.get('application.log_level', 'INFO')
    
    def get_theme(self) -> str:
        """テーマの取得"""
        return self.get('gui.theme', 'dark')
    
    def get_language(self) -> str:
        """言語設定の取得"""
        return self.get('gui.language', 'ja')