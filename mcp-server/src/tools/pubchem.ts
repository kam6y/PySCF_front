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
          text: `**Compound Search Results**

**Compound Information:**
- **Name:** ${compound.iupac_name || 'N/A'}
- **Molecular Formula:** ${compound.molecular_formula || 'N/A'}
- **Molecular Weight:** ${formatMolecularWeight(compound.molecular_weight)} g/mol
- **PubChem CID:** ${compound.cid}
- **Number of Atoms:** ${atomCount}

**Synonyms:** ${compound.synonyms.slice(0, 5).join(', ')}${compound.synonyms.length > 5 ? '...' : ''}

**XYZ Structure Data:**
\`\`\`
${response.data.xyz}
\`\`\`

This 3D structure data can be used for quantum chemistry calculations.`,
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
**Debug Information:**
- HTTP Status: ${axiosError.response?.status || 'N/A'}
- Response Data: ${JSON.stringify(axiosError.response?.data, null, 2) || 'N/A'}`;
    }
    
    return {
      content: [
        {
          type: 'text',
          text: `❌ PubChem search error occurred: ${errorMessage}

**Input Parameters:**
- Query: "${args.query}"
- Search Type: "${args.search_type || 'name'}"${debugInfo}

**Possible Causes:**
- Misspelled compound name
- Compound not registered in PubChem
- Network connection issues
- Server error

**Solutions:**
- Check the spelling of the compound name
- Try synonyms or IUPAC names
- If you know the CID, set search_type to 'cid'`,
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
          text: `**SMILES Conversion Results**

**Input SMILES:** \`${args.smiles}\`

**Generated 3D XYZ Structure:**
\`\`\`
${response.data.xyz}
\`\`\`

This 3D structure data can be used for quantum chemistry calculations. Performing geometry optimization will provide more accurate molecular geometry.`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ SMILES conversion error occurred: ${errorMessage}

**Possible Causes:**
- Invalid SMILES string
- Unsupported molecular structure
- Failed 3D coordinate generation

**Solutions:**
- Check the syntax of the SMILES string
- Try with simpler molecular structures first
- If stereochemistry is not explicit, try different SMILES notation

**SMILES Examples:**
- Water: O
- Methane: C
- Benzene: c1ccccc1
- Ethanol: CCO`,
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
            text: `❌ **Invalid XYZ format**

**Error:** ${data.error || 'Unknown error'}

**Correct XYZ format:**
\`\`\`
[Number of atoms]
[Title line (optional)]
[Element symbol] [X coordinate] [Y coordinate] [Z coordinate]
[Element symbol] [X coordinate] [Y coordinate] [Z coordinate]
...
\`\`\`

**Example (Water molecule):**
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
          text: `✅ **XYZ format is valid**

**Molecular Information:**
- **Number of atoms:** ${data.num_atoms}
- **Molecular formula:** ${formulaText}
- **Title:** ${data.title || '(none)'}

**Atomic composition:**
${Object.entries(elementCounts)
  .map(([element, count]) => `- ${element}: ${count} atoms`)
  .join('\n')}

**Atomic coordinates:**
\`\`\`
${data.atoms?.map((atom, i) =>
  `${(i + 1).toString().padStart(2)} ${atom.element.padEnd(2)} ${atom.x.toFixed(6).padStart(10)} ${atom.y.toFixed(6).padStart(10)} ${atom.z.toFixed(6).padStart(10)}`
).join('\n') || ''}
\`\`\`

This molecular structure can be used for quantum chemistry calculations.`,
        },
      ],
    };
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
    return {
      content: [
        {
          type: 'text',
          text: `❌ XYZ validation error occurred: ${errorMessage}

A server error occurred during XYZ format validation. Please ensure the server is running.`,
        },
      ],
      isError: true,
    };
  }
}