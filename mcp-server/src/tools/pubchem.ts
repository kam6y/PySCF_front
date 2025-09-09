import { Tool } from '@modelcontextprotocol/sdk/types.js';
import { PySCFApiClient, PubChemSearchRequest, SMILESConvertRequest, XYZValidateRequest } from '../client.js';

export const searchPubChemTool: Tool = {
  name: 'searchPubChem',
  description: 'Search PubChem database for chemical compounds and retrieve their 3D molecular structure in XYZ format',
  inputSchema: {
    type: 'object',
    properties: {
      query: {
        type: 'string',
        description: 'Search query for PubChem (compound name, CID, or molecular formula)',
        minLength: 1,
      },
      search_type: {
        type: 'string',
        enum: ['name', 'cid', 'formula'],
        description: 'Type of search to perform',
        default: 'name',
      },
    },
    required: ['query'],
  },
};

export async function handleSearchPubChem(
  args: { query: string; search_type?: 'name' | 'cid' | 'formula' },
  client: PySCFApiClient
) {
  try {
    const request: PubChemSearchRequest = {
      query: args.query,
      search_type: args.search_type || 'name',
    };

    const response = await client.searchPubChem(request);

    if (!response.success) {
      throw new Error(`PubChem search failed`);
    }

    const compound = response.data.compound_info;
    const atomCount = response.data.atom_count;

    // Safely format molecular weight with type checking
    const formatMolecularWeight = (weight: any): string => {
      if (typeof weight === 'number' && !isNaN(weight)) {
        return weight.toFixed(2);
      }
      if (typeof weight === 'string' && !isNaN(parseFloat(weight))) {
        return parseFloat(weight).toFixed(2);
      }
      return 'N/A';
    };

    return {
      content: [
        {
          type: 'text',
          text: `**化合物検索結果**

**化合物情報:**
- **名前:** ${compound.iupac_name || 'N/A'}
- **分子式:** ${compound.molecular_formula || 'N/A'}
- **分子量:** ${formatMolecularWeight(compound.molecular_weight)} g/mol
- **PubChem CID:** ${compound.cid}
- **原子数:** ${atomCount}

**同義語:** ${compound.synonyms.slice(0, 5).join(', ')}${compound.synonyms.length > 5 ? '...' : ''}

**XYZ構造データ:**
\`\`\`
${response.data.xyz}
\`\`\`

この3D構造データは量子化学計算で使用できます。`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    
    // Enhanced error details for debugging
    let debugInfo = '';
    if (error instanceof Error && 'response' in error) {
      const axiosError = error as any;
      debugInfo = `
**デバッグ情報:**
- HTTPステータス: ${axiosError.response?.status || 'N/A'}
- レスポンスデータ: ${JSON.stringify(axiosError.response?.data, null, 2) || 'N/A'}`;
    }
    
    return {
      content: [
        {
          type: 'text',
          text: `❌ PubChem検索でエラーが発生しました: ${errorMessage}

**入力パラメータ:**
- クエリ: "${args.query}"
- 検索タイプ: "${args.search_type || 'name'}"${debugInfo}

**可能な原因:**
- 化合物名の入力ミス
- PubChemに登録されていない化合物
- ネットワーク接続の問題
- サーバーエラー

**解決方法:**
- 化合物名のスペルを確認してください
- 同義語やIUPAC名を試してください  
- CIDが分かる場合はsearch_typeを'cid'に設定してください`,
        },
      ],
      isError: true,
    };
  }
}

export const convertSmilesTool: Tool = {
  name: 'convertSmiles',
  description: 'Convert SMILES string to 3D molecular structure in XYZ format for quantum chemistry calculations',
  inputSchema: {
    type: 'object',
    properties: {
      smiles: {
        type: 'string',
        description: 'SMILES string representing the molecular structure',
        minLength: 1,
      },
    },
    required: ['smiles'],
  },
};

export async function handleConvertSmiles(args: { smiles: string }, client: PySCFApiClient) {
  try {
    const request: SMILESConvertRequest = {
      smiles: args.smiles,
    };

    const response = await client.convertSmiles(request);

    if (!response.success) {
      throw new Error(`SMILES conversion failed`);
    }

    return {
      content: [
        {
          type: 'text',
          text: `**SMILES変換結果**

**入力SMILES:** \`${args.smiles}\`

**生成された3D XYZ構造:**
\`\`\`
${response.data.xyz}
\`\`\`

この3D構造データは量子化学計算で使用できます。構造最適化を行うとより正確な分子形状が得られます。`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ SMILES変換でエラーが発生しました: ${errorMessage}

**可能な原因:**
- 無効なSMILES文字列
- サポートされていない分子構造
- 3D座標生成の失敗

**解決方法:**
- SMILES文字列の構文を確認してください
- より単純な分子構造から試してください
- 立体化学が明示的でない場合は、異なるSMILES表記を試してください

**SMILES例:**
- 水: O
- メタン: C
- ベンゼン: c1ccccc1
- エタノール: CCO`,
        },
      ],
      isError: true,
    };
  }
}

export const validateXyzTool: Tool = {
  name: 'validateXYZ',
  description: 'Validate XYZ format molecular structure data and extract molecular information',
  inputSchema: {
    type: 'object',
    properties: {
      xyz: {
        type: 'string',
        description: 'XYZ format molecular structure data to validate',
        minLength: 1,
      },
    },
    required: ['xyz'],
  },
};

export async function handleValidateXYZ(args: { xyz: string }, client: PySCFApiClient) {
  try {
    const request: XYZValidateRequest = {
      xyz: args.xyz,
    };

    const response = await client.validateXYZ(request);

    if (!response.success) {
      throw new Error(`XYZ validation failed`);
    }

    const data = response.data;

    if (!data.valid) {
      return {
        content: [
          {
            type: 'text',
            text: `❌ **XYZ形式が無効です**

**エラー:** ${data.error || '不明なエラー'}

**正しいXYZ形式:**
\`\`\`
[原子数]
[タイトル行（省略可能）]
[元素記号] [X座標] [Y座標] [Z座標]
[元素記号] [X座標] [Y座標] [Z座標]
...
\`\`\`

**例（水分子）:**
\`\`\`
3
Water molecule
O  0.000000  0.000000  0.000000
H  0.757000  0.587000  0.000000
H -0.757000  0.587000  0.000000
\`\`\``,
          },
        ],
        isError: true,
      };
    }

    // Group atoms by element
    const elementCounts: { [element: string]: number } = {};
    data.atoms?.forEach(atom => {
      elementCounts[atom.element] = (elementCounts[atom.element] || 0) + 1;
    });

    const formulaText = Object.entries(elementCounts)
      .sort()
      .map(([element, count]) => count > 1 ? `${element}${count}` : element)
      .join('');

    return {
      content: [
        {
          type: 'text',
          text: `✅ **XYZ形式は有効です**

**分子情報:**
- **原子数:** ${data.num_atoms}
- **分子式:** ${formulaText}
- **タイトル:** ${data.title || '(なし)'}

**原子構成:**
${Object.entries(elementCounts)
  .map(([element, count]) => `- ${element}: ${count}個`)
  .join('\n')}

**原子座標:**
\`\`\`
${data.atoms?.map((atom, i) => 
  `${(i + 1).toString().padStart(2)} ${atom.element.padEnd(2)} ${atom.x.toFixed(6).padStart(10)} ${atom.y.toFixed(6).padStart(10)} ${atom.z.toFixed(6).padStart(10)}`
).join('\n') || ''}
\`\`\`

この分子構造は量子化学計算で使用できます。`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ XYZ検証でエラーが発生しました: ${errorMessage}

XYZ形式の検証中にサーバーエラーが発生しました。サーバーが起動していることを確認してください。`,
        },
      ],
      isError: true,
    };
  }
}