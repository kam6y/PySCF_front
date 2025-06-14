# MCPサーバー設定（オプション機能）

## 概要
MCPサーバーはPySCF_Frontのオプション機能で、Claude Desktopなどの外部アプリケーションから自然言語でPySCF_Frontを操作できるようにします。

## MCPサーバーの実装

### backend/mcp_server/server.py
```python
"""
PySCF_Front MCPサーバー実装
Claude Desktopなどから接続可能なMCPサーバー
"""
import asyncio
import json
import logging
from typing import Dict, Any, List, Optional

from mcp import Server, Tool, Resource
from mcp.server.stdio import stdio_server
from fastmcp import FastMCP

from core.calculation_engine import CalculationEngine
from core.molecule_manager import MoleculeManager
from database.repository import Repository
from utils.logger import setup_logger

logger = setup_logger(__name__)


class PySCFMCPServer:
    """PySCF_Front MCPサーバー"""
    
    def __init__(self, settings):
        self.settings = settings
        self.mcp = FastMCP("PySCF_Front MCP Server")
        self.calc_engine = CalculationEngine(settings)
        self.mol_manager = MoleculeManager()
        self.repository = Repository(settings.database_url)
        
        # ツールの登録
        self._register_tools()
    
    def _register_tools(self):
        """MCPツールの登録"""
        
        @self.mcp.tool()
        async def create_molecule_from_text(description: str) -> Dict[str, Any]:
            """
            自然言語の説明から分子構造を生成
            
            Args:
                description: 分子の説明（例: "水分子", "ベンゼン環", "メタノール"）
            
            Returns:
                生成された分子情報
            """
            try:
                # 分子名から構造を生成
                molecule = await self.mol_manager.create_from_description(description)
                
                # データベースに保存
                mol_id = await self.repository.save_molecule(molecule)
                
                return {
                    "success": True,
                    "molecule_id": mol_id,
                    "name": molecule.name,
                    "formula": molecule.formula,
                    "atoms": molecule.to_xyz(),
                    "charge": molecule.charge,
                    "multiplicity": molecule.multiplicity
                }
            except Exception as e:
                logger.error(f"Failed to create molecule: {e}")
                return {
                    "success": False,
                    "error": str(e)
                }
        
        @self.mcp.tool()
        async def get_calculation_recommendation(
            molecule_id: str,
            purpose: str = "energy"
        ) -> Dict[str, Any]:
            """
            分子と計算目的に基づいて最適な計算手法を推奨
            
            Args:
                molecule_id: 分子ID
                purpose: 計算目的 (energy, optimization, frequency, excited_state)
            
            Returns:
                推奨される計算設定
            """
            try:
                # 分子情報を取得
                molecule = await self.repository.get_molecule(molecule_id)
                if not molecule:
                    return {"success": False, "error": "Molecule not found"}
                
                atom_count = len(molecule.atoms)
                
                # 推奨設定の決定
                if purpose == "energy":
                    if atom_count < 50:
                        method = "B3LYP"
                        basis = "6-311G(d,p)"
                        reason = "中規模分子のエネルギー計算に適したDFT汎関数と基底関数"
                    else:
                        method = "PBE"
                        basis = "def2-SVP"
                        reason = "大規模分子に適した効率的なDFT計算"
                
                elif purpose == "optimization":
                    if atom_count < 30:
                        method = "B3LYP"
                        basis = "6-31G(d)"
                        reason = "構造最適化に適したバランスの良い設定"
                    else:
                        method = "PBE"
                        basis = "def2-SVP"
                        reason = "大規模系の構造最適化に適した設定"
                
                elif purpose == "frequency":
                    method = "B3LYP"
                    basis = "6-31G(d)"
                    reason = "振動解析に必要な精度と計算効率のバランス"
                
                elif purpose == "excited_state":
                    method = "TD-B3LYP"
                    basis = "6-31+G(d)"
                    reason = "励起状態計算に適したTDDFT設定"
                
                else:
                    method = "HF"
                    basis = "6-31G"
                    reason = "基本的な計算設定"
                
                return {
                    "success": True,
                    "method": method,
                    "basis_set": basis,
                    "reason": reason,
                    "estimated_time": self._estimate_calculation_time(atom_count, method)
                }
                
            except Exception as e:
                logger.error(f"Failed to get recommendation: {e}")
                return {
                    "success": False,
                    "error": str(e)
                }
        
        @self.mcp.tool()
        async def start_calculation(
            molecule_id: str,
            method: Optional[str] = None,
            basis_set: Optional[str] = None,
            **kwargs
        ) -> Dict[str, Any]:
            """
            計算を開始
            
            Args:
                molecule_id: 分子ID
                method: 計算手法（省略時は推奨設定を使用）
                basis_set: 基底関数（省略時は推奨設定を使用）
                **kwargs: その他の計算パラメータ
            
            Returns:
                計算ジョブ情報
            """
            try:
                # 推奨設定を取得（指定がない場合）
                if not method or not basis_set:
                    recommendation = await self.get_calculation_recommendation(
                        molecule_id, kwargs.get("purpose", "energy")
                    )
                    method = method or recommendation["method"]
                    basis_set = basis_set or recommendation["basis_set"]
                
                # 計算ジョブを作成
                job_id = await self.calc_engine.submit_calculation(
                    molecule_id=molecule_id,
                    method=method,
                    basis_set=basis_set,
                    parameters=kwargs
                )
                
                return {
                    "success": True,
                    "job_id": job_id,
                    "status": "submitted",
                    "method": method,
                    "basis_set": basis_set
                }
                
            except Exception as e:
                logger.error(f"Failed to start calculation: {e}")
                return {
                    "success": False,
                    "error": str(e)
                }
        
        @self.mcp.tool()
        async def get_job_status(job_id: str) -> Dict[str, Any]:
            """
            計算ジョブの状態を確認
            
            Args:
                job_id: ジョブID
            
            Returns:
                ジョブの状態情報
            """
            try:
                status = await self.calc_engine.get_job_status(job_id)
                return {
                    "success": True,
                    "job_id": job_id,
                    "status": status.status,
                    "progress": status.progress,
                    "message": status.message,
                    "started_at": status.started_at,
                    "updated_at": status.updated_at
                }
            except Exception as e:
                logger.error(f"Failed to get job status: {e}")
                return {
                    "success": False,
                    "error": str(e)
                }
        
        @self.mcp.tool()
        async def get_calculation_results(
            job_id: str,
            include_details: bool = True
        ) -> Dict[str, Any]:
            """
            計算結果を取得
            
            Args:
                job_id: ジョブID
                include_details: 詳細情報を含めるか
            
            Returns:
                計算結果
            """
            try:
                results = await self.repository.get_results(job_id)
                if not results:
                    return {"success": False, "error": "Results not found"}
                
                response = {
                    "success": True,
                    "job_id": job_id,
                    "total_energy": results.total_energy,
                    "energy_unit": "Hartree",
                    "converged": results.converged,
                    "calculation_time": results.calculation_time
                }
                
                if include_details:
                    response.update({
                        "homo_energy": results.homo_energy,
                        "lumo_energy": results.lumo_energy,
                        "homo_lumo_gap": results.homo_lumo_gap,
                        "dipole_moment": results.dipole_moment,
                        "orbital_energies": results.orbital_energies[:10]  # 最初の10個
                    })
                
                return response
                
            except Exception as e:
                logger.error(f"Failed to get results: {e}")
                return {
                    "success": False,
                    "error": str(e)
                }
        
        @self.mcp.tool()
        async def interpret_results(
            job_id: str,
            level: str = "basic"
        ) -> Dict[str, Any]:
            """
            計算結果を解釈して説明
            
            Args:
                job_id: ジョブID
                level: 解釈レベル (basic, detailed, educational)
            
            Returns:
                結果の解釈
            """
            try:
                results = await self.repository.get_results(job_id)
                if not results:
                    return {"success": False, "error": "Results not found"}
                
                interpretation = {
                    "success": True,
                    "job_id": job_id,
                    "summary": ""
                }
                
                # 基本的な解釈
                if results.converged:
                    interpretation["summary"] = "計算は正常に収束しました。"
                else:
                    interpretation["summary"] = "計算は収束しませんでした。設定の見直しが必要です。"
                
                # エネルギーの解釈
                energy_ev = results.total_energy * 27.2114  # Hartree to eV
                interpretation["energy_interpretation"] = (
                    f"全エネルギーは {results.total_energy:.6f} Hartree "
                    f"({energy_ev:.2f} eV) です。"
                )
                
                # HOMO-LUMOギャップの解釈
                if results.homo_lumo_gap:
                    gap_ev = results.homo_lumo_gap * 27.2114
                    interpretation["homo_lumo_interpretation"] = (
                        f"HOMO-LUMOギャップは {gap_ev:.2f} eV です。"
                    )
                    
                    if gap_ev < 1.5:
                        interpretation["homo_lumo_interpretation"] += (
                            " これは比較的小さなギャップで、導電性や反応性が高い可能性があります。"
                        )
                    elif gap_ev > 5.0:
                        interpretation["homo_lumo_interpretation"] += (
                            " これは大きなギャップで、化学的に安定な分子であることを示しています。"
                        )
                    else:
                        interpretation["homo_lumo_interpretation"] += (
                            " これは典型的な有機分子のギャップ値です。"
                        )
                
                # 詳細レベルの解釈
                if level in ["detailed", "educational"]:
                    interpretation["detailed_insights"] = []
                    
                    # 双極子モーメントの解釈
                    if results.dipole_moment:
                        dipole_debye = results.dipole_moment * 2.5417  # a.u. to Debye
                        insight = f"双極子モーメントは {dipole_debye:.2f} Debye です。"
                        if dipole_debye > 2.0:
                            insight += " これは極性が高い分子であることを示しています。"
                        interpretation["detailed_insights"].append(insight)
                
                # 教育的レベルの解釈
                if level == "educational":
                    interpretation["educational_notes"] = [
                        "Hartreeは原子単位系でのエネルギー単位です。",
                        "HOMO（最高被占軌道）とLUMO（最低空軌道）のエネルギー差は、分子の化学的性質を理解する上で重要です。",
                        "双極子モーメントは分子内の電荷分布の偏りを表し、分子の極性を示します。"
                    ]
                
                return interpretation
                
            except Exception as e:
                logger.error(f"Failed to interpret results: {e}")
                return {
                    "success": False,
                    "error": str(e)
                }
    
    def _estimate_calculation_time(self, atom_count: int, method: str) -> str:
        """計算時間の推定"""
        # 簡単な推定ロジック
        base_time = 1.0  # 分
        
        if method.startswith("TD-"):
            base_time *= 3
        elif method in ["MP2", "CCSD"]:
            base_time *= 5
        elif method == "CCSD(T)":
            base_time *= 10
        
        time_minutes = base_time * (atom_count / 10) ** 2
        
        if time_minutes < 60:
            return f"{int(time_minutes)}分"
        else:
            return f"{time_minutes / 60:.1f}時間"
    
    async def start(self):
        """MCPサーバーの起動"""
        logger.info(f"Starting MCP server on port {self.settings.mcp_port}")
        
        # stdioサーバーとして起動（Claude Desktopとの通信用）
        async with stdio_server() as server:
            await server.run(
                self.mcp.get_handler(),
                server_name="PySCF_Front MCP Server",
                server_version="1.0.0"
            )


async def start_mcp_server(settings):
    """MCPサーバーを起動"""
    server = PySCFMCPServer(settings)
    await server.start()
```

### backend/mcp_server/tools.py
```python
"""
MCPツールの定義と実装
"""
from typing import List, Dict, Any
from mcp import Tool


def get_pyscf_tools() -> List[Tool]:
    """PySCF_Front用のMCPツール定義を返す"""
    
    tools = [
        Tool(
            name="create_molecule_from_text",
            description="自然言語の説明から分子構造を生成します",
            input_schema={
                "type": "object",
                "properties": {
                    "description": {
                        "type": "string",
                        "description": "分子の説明（例：水分子、ベンゼン環、エタノール）"
                    }
                },
                "required": ["description"]
            }
        ),
        
        Tool(
            name="get_calculation_recommendation",
            description="分子と計算目的に基づいて最適な計算手法を推奨します",
            input_schema={
                "type": "object",
                "properties": {
                    "molecule_id": {
                        "type": "string",
                        "description": "分子ID"
                    },
                    "purpose": {
                        "type": "string",
                        "enum": ["energy", "optimization", "frequency", "excited_state"],
                        "description": "計算目的",
                        "default": "energy"
                    }
                },
                "required": ["molecule_id"]
            }
        ),
        
        Tool(
            name="start_calculation",
            description="量子化学計算を開始します",
            input_schema={
                "type": "object",
                "properties": {
                    "molecule_id": {
                        "type": "string",
                        "description": "分子ID"
                    },
                    "method": {
                        "type": "string",
                        "description": "計算手法（省略時は推奨設定を使用）"
                    },
                    "basis_set": {
                        "type": "string",
                        "description": "基底関数（省略時は推奨設定を使用）"
                    }
                },
                "required": ["molecule_id"]
            }
        ),
        
        Tool(
            name="get_job_status",
            description="計算ジョブの現在の状態を確認します",
            input_schema={
                "type": "object",
                "properties": {
                    "job_id": {
                        "type": "string",
                        "description": "ジョブID"
                    }
                },
                "required": ["job_id"]
            }
        ),
        
        Tool(
            name="get_calculation_results",
            description="計算結果を取得します",
            input_schema={
                "type": "object",
                "properties": {
                    "job_id": {
                        "type": "string",
                        "description": "ジョブID"
                    },
                    "include_details": {
                        "type": "boolean",
                        "description": "詳細情報を含めるか",
                        "default": True
                    }
                },
                "required": ["job_id"]
            }
        ),
        
        Tool(
            name="interpret_results",
            description="計算結果を解釈して分かりやすく説明します",
            input_schema={
                "type": "object",
                "properties": {
                    "job_id": {
                        "type": "string",
                        "description": "ジョブID"
                    },
                    "level": {
                        "type": "string",
                        "enum": ["basic", "detailed", "educational"],
                        "description": "解釈の詳細レベル",
                        "default": "basic"
                    }
                },
                "required": ["job_id"]
            }
        )
    ]
    
    return tools
```

### backend/config/mcp_config.yaml
```yaml
# MCP サーバー設定ファイル
mcp_server:
  # MCPサーバーの有効/無効
  enabled: false  # デフォルトは無効
  
  # 接続設定
  host: "localhost"
  port: 50053
  
  # セキュリティ設定
  security:
    # リモート接続を許可するか
    allow_remote: false
    
    # 認証設定
    authentication:
      enabled: false
      method: "token"  # token, oauth2
      token_file: "/etc/pyscf_front/mcp_token"
    
    # 接続制限
    rate_limiting:
      enabled: true
      max_requests_per_minute: 60
      max_concurrent_requests: 10
    
    # 許可するIPアドレス
    allowed_ips:
      - "127.0.0.1"
      - "::1"
  
  # 利用可能な機能
  allowed_tools:
    - "create_molecule_from_text"
    - "get_calculation_recommendation"
    - "start_calculation"
    - "get_job_status"
    - "get_calculation_results"
    - "interpret_results"
  
  # リソース制限
  resource_limits:
    max_molecule_size: 1000  # 最大原子数
    max_calculation_time: 3600  # 秒
    max_memory_per_job: 8192  # MB
  
  # ログ設定
  logging:
    level: "INFO"
    file: "/var/log/pyscf_front/mcp_server.log"
```

## Claude Desktop設定

### 1. Claude Desktop設定ファイル
Claude Desktopの設定ファイル（通常は `~/Library/Application Support/Claude/claude_desktop_config.json` on macOS）に以下を追加：

```json
{
  "mcpServers": {
    "pyscf-front": {
      "command": "docker",
      "args": [
        "exec", 
        "-i", 
        "pyscf_backend", 
        "python", 
        "-m", 
        "mcp_server.server"
      ],
      "env": {
        "PYTHONPATH": "/app",
        "MCP_SERVER_MODE": "stdio"
      }
    }
  }
}
```

### 2. ローカル開発用の設定
開発時は以下の設定も使用可能：

```json
{
  "mcpServers": {
    "pyscf-front-dev": {
      "command": "python",
      "args": [
        "-m",
        "mcp_server.server"
      ],
      "cwd": "/path/to/pyscf_front/backend",
      "env": {
        "PYTHONPATH": "/path/to/pyscf_front/backend",
        "MCP_SERVER_ENABLED": "true",
        "MYSQL_PASSWORD": "your_password"
      }
    }
  }
}
```

## MCPサーバーの起動方法

### 1. Docker Composeで起動
```bash
# MCP有効化オプション付きで起動
docker-compose -f docker-compose.yml -f docker-compose.mcp.yml up -d
```

### 2. 環境変数で有効化
```bash
# .envファイルに追加
MCP_SERVER_ENABLED=true

# 通常のDocker Compose起動
docker-compose up -d
```

### 3. 手動起動（開発時）
```bash
cd backend
export MCP_SERVER_ENABLED=true
python -m mcp_server.server
```

## 使用例

### Claude Desktopでの対話例

```
ユーザー: PySCF_Frontを使って水分子の計算をしたい

Claude: 水分子の量子化学計算を行います。まず分子構造を作成します。

[create_molecule_from_text("水分子")]

水分子（H2O）の構造を作成しました。
- 分子ID: mol_12345
- 化学式: H2O
- 原子数: 3
- 電荷: 0
- スピン多重度: 1

次に、この分子に適した計算条件を確認します。

[get_calculation_recommendation("mol_12345", "energy")]

エネルギー計算には以下の設定を推奨します：
- 計算手法: B3LYP
- 基底関数: 6-311G(d,p)
- 推定計算時間: 約1分

この設定で計算を開始してよろしいですか？

ユーザー: はい、お願いします

Claude: 計算を開始します。

[start_calculation("mol_12345")]

計算が開始されました。
- ジョブID: job_67890
- 状態: 実行中

進捗を確認します...

[get_job_status("job_67890")]

計算が完了しました！結果を取得します。

[get_calculation_results("job_67890")]

計算結果：
- 全エネルギー: -76.4183 Hartree
- HOMO: -0.3912 Hartree
- LUMO: 0.0523 Hartree
- HOMO-LUMOギャップ: 0.4435 Hartree (12.07 eV)
- 双極子モーメント: 0.8162 a.u. (2.08 Debye)

[interpret_results("job_67890", "detailed")]

結果の解釈：
- 計算は正常に収束しました
- HOMO-LUMOギャップが12.07 eVと大きく、水分子は化学的に安定です
- 双極子モーメントが2.08 Debyeと比較的大きく、極性分子であることを示しています
- これらの値は文献値とよく一致しており、計算が正しく実行されたことを確認できます
```

## セキュリティに関する注意事項

1. **デフォルトでは無効**
   - MCPサーバーはオプション機能で、明示的に有効化する必要があります

2. **ローカル接続のみ**
   - デフォルトではlocalhostからの接続のみ許可

3. **認証オプション**
   - 必要に応じてトークン認証を有効化できます

4. **リソース制限**
   - 計算時間やメモリ使用量に制限を設定可能

5. **操作制限**
   - 利用可能なツールを制限できます

これらの設定により、PySCF_FrontをClaude Desktopから安全に操作できます。