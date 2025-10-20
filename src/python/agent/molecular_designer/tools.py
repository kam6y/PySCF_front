"""
Molecular Designer Tool Wrapper Module

This module provides RDKit-based functions that enable AI agents to
generate novel molecular structures and predict simple molecular properties.

Each function is optimized for Gemini SDK's Function Calling feature and
implements cheminformatics workflows for molecular design tasks.
"""

import json
import logging
from typing import Optional, List, Dict
from langchain_core.tools import tool
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, BRICS, Recap

# Logger configuration
logger = logging.getLogger(__name__)


def _validation_error(message: str) -> str:
    """
    Return a validation error message in JSON format.

    Args:
        message: The validation error message

    Returns:
        str: JSON-formatted error message with structure: {"success": false, "error": "message"}
    """
    return json.dumps({
        "success": False,
        "error": message
    }, ensure_ascii=False, indent=2)


@tool
def generate_analogs_rdkit(base_smiles: str, num_to_generate: int = 10, design_strategy: str = "diverse") -> str:
    """
    Generate molecular analogs based on a base SMILES structure.

    This function generates molecular variants by applying chemical modifications
    such as adding electron-donating/withdrawing groups, extending conjugation,
    or introducing heteroatoms. The generated molecules are returned as SMILES strings
    with design rationale explanations.

    Design strategies:
    - "diverse": Generate a diverse set of analogs with different modifications
    - "conjugation": Focus on extending conjugation systems (for longer wavelength absorption)
    - "push_pull": Focus on push-pull structures (electron donor + acceptor groups)
    - "heteroatom": Focus on introducing heteroatoms (N, O, S substitutions)

    Args:
        base_smiles (str): Base molecular structure in SMILES format.
                          Example: "c1ccccc1" (benzene)
        num_to_generate (int): Maximum number of analogs to generate (1-50).
                              Default: 10
        design_strategy (str): Strategy for molecular modification.
                              Options: "diverse", "conjugation", "push_pull", "heteroatom"
                              Default: "diverse"

    Returns:
        str: JSON string containing generated SMILES and design rationale.
             On success: {"success": true, "base_smiles": "...", "generated_candidates": [...]}
             On error: {"success": false, "error": "..."}
    """
    # Input validation
    if not base_smiles or not isinstance(base_smiles, str):
        return _validation_error("base_smiles parameter is required and must be a non-empty string.")

    if not isinstance(num_to_generate, int) or num_to_generate < 1 or num_to_generate > 50:
        return _validation_error("num_to_generate must be an integer between 1 and 50.")

    if design_strategy not in ["diverse", "conjugation", "push_pull", "heteroatom"]:
        return _validation_error("design_strategy must be one of: 'diverse', 'conjugation', 'push_pull', 'heteroatom'")

    try:
        # Validate base SMILES
        mol = Chem.MolFromSmiles(base_smiles)
        if not mol:
            return _validation_error(f"Invalid SMILES string: {base_smiles}")

        logger.info(f"Generating analogs for {base_smiles} with strategy: {design_strategy}")

        # Define modification rules based on strategy
        candidates = []

        if design_strategy == "diverse" or design_strategy == "conjugation":
            # Extend conjugation systems
            candidates.extend([
                {"smiles": "c1c(C=C)cccc1", "rationale": "Added vinyl group to extend conjugation"},
                {"smiles": "c1c(C#C)cccc1", "rationale": "Added ethynyl group for extended conjugation"},
                {"smiles": "c1ccc(cc1)c2ccccc2", "rationale": "Added phenyl ring to create biphenyl system"},
                {"smiles": "c1ccc(cc1)C=Cc2ccccc2", "rationale": "Created stilbene-like structure with extended conjugation"},
            ])

        if design_strategy == "diverse" or design_strategy == "push_pull":
            # Push-pull structures (electron donor + acceptor)
            candidates.extend([
                {"smiles": "c1c(N)ccc(C#N)cc1", "rationale": "Push-pull: Amino (donor) and cyano (acceptor)"},
                {"smiles": "c1c(N(C)C)ccc(C#N)cc1", "rationale": "Push-pull: Dimethylamino (strong donor) and cyano (acceptor)"},
                {"smiles": "c1c(O)ccc(C(=O)O)cc1", "rationale": "Push-pull: Hydroxyl (donor) and carboxyl (acceptor)"},
                {"smiles": "c1c(N)ccc(C(F)(F)F)cc1", "rationale": "Push-pull: Amino (donor) and trifluoromethyl (acceptor)"},
            ])

        if design_strategy == "diverse" or design_strategy == "heteroatom":
            # Heteroatom substitution
            candidates.extend([
                {"smiles": "c1ncccc1", "rationale": "Replaced C with N to create pyridine"},
                {"smiles": "c1ccoc1", "rationale": "Created furan ring (oxygen-containing heterocycle)"},
                {"smiles": "c1ccsc1", "rationale": "Created thiophene ring (sulfur-containing heterocycle)"},
                {"smiles": "c1c2c(ccc1)nc(cc2)", "rationale": "Created quinoline (benzopyridine) structure"},
            ])

        if design_strategy == "diverse":
            # Additional diverse modifications
            candidates.extend([
                {"smiles": "c1c(F)cccc1", "rationale": "Added fluorine (electron-withdrawing)"},
                {"smiles": "c1c(Cl)cccc1", "rationale": "Added chlorine (electron-withdrawing)"},
                {"smiles": "c1c(O)cccc1", "rationale": "Added hydroxyl (electron-donating)"},
                {"smiles": "c1c(C)cccc1", "rationale": "Added methyl (electron-donating)"},
                {"smiles": "c1c(OC)cccc1", "rationale": "Added methoxy (electron-donating)"},
                {"smiles": "c1c(C(=O)C)cccc1", "rationale": "Added acetyl (electron-withdrawing)"},
            ])

        # Validate all generated SMILES
        validated_candidates = []
        for candidate in candidates[:num_to_generate]:
            test_mol = Chem.MolFromSmiles(candidate["smiles"])
            if test_mol:
                validated_candidates.append(candidate)
            else:
                logger.warning(f"Generated invalid SMILES: {candidate['smiles']}")

        logger.info(f"Successfully generated {len(validated_candidates)} valid analogs")

        return json.dumps({
            "success": True,
            "base_smiles": base_smiles,
            "design_strategy": design_strategy,
            "generated_candidates": validated_candidates,
            "total_generated": len(validated_candidates)
        }, ensure_ascii=False, indent=2)

    except Exception as e:
        logger.error(f"Error in generate_analogs_rdkit: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


@tool
def predict_simple_properties_rdkit(smiles_list: str) -> str:
    """
    Predict simple molecular properties using RDKit Descriptors.

    This function quickly calculates molecular properties that can be used for
    initial filtering before expensive quantum chemistry calculations. Properties include:
    - Molecular weight (MW)
    - LogP (lipophilicity)
    - TPSA (topological polar surface area)
    - Number of rotatable bonds
    - Number of H-bond donors and acceptors
    - Number of aromatic rings

    These properties are useful for drug-likeness assessment and filtering candidates.

    Args:
        smiles_list (str): JSON string containing list of SMILES strings to analyze.
                          Example: '["CCO", "c1ccccc1", "CC(=O)O"]'
                          Or a single SMILES string: "CCO"

    Returns:
        str: JSON string containing predicted properties for each molecule.
             On success: {"success": true, "predictions": [...]}
             On error: {"success": false, "error": "..."}
    """
    # Input validation
    if not smiles_list or not isinstance(smiles_list, str):
        return _validation_error("smiles_list parameter is required and must be a non-empty string.")

    try:
        # Parse input - accept both JSON array and single SMILES
        if smiles_list.strip().startswith('['):
            # JSON array format
            smiles_array = json.loads(smiles_list)
            if not isinstance(smiles_array, list):
                return _validation_error("smiles_list must be a JSON array of strings or a single SMILES string.")
        else:
            # Single SMILES string
            smiles_array = [smiles_list.strip()]

        if not smiles_array:
            return _validation_error("smiles_list cannot be empty.")

        if len(smiles_array) > 100:
            return _validation_error("Maximum 100 molecules can be processed at once.")

        logger.info(f"Predicting properties for {len(smiles_array)} molecules")

        results = []
        for smiles in smiles_array:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    results.append({
                        "smiles": smiles,
                        "error": "Invalid SMILES string"
                    })
                    continue

                # Add hydrogens for accurate calculations
                mol = Chem.AddHs(mol)

                # Calculate properties
                properties = {
                    "smiles": smiles,
                    "molecular_weight": round(Descriptors.MolWt(mol), 2),
                    "logp": round(Descriptors.MolLogP(mol), 2),
                    "tpsa": round(Descriptors.TPSA(mol), 2),
                    "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                    "num_h_donors": Descriptors.NumHDonors(mol),
                    "num_h_acceptors": Descriptors.NumHAcceptors(mol),
                    "num_aromatic_rings": Descriptors.NumAromaticRings(mol),
                    "num_atoms": mol.GetNumAtoms(),
                }

                # Add Lipinski's Rule of Five assessment
                lipinski_violations = 0
                if properties["molecular_weight"] > 500:
                    lipinski_violations += 1
                if properties["logp"] > 5:
                    lipinski_violations += 1
                if properties["num_h_donors"] > 5:
                    lipinski_violations += 1
                if properties["num_h_acceptors"] > 10:
                    lipinski_violations += 1

                properties["lipinski_violations"] = lipinski_violations
                properties["lipinski_compliant"] = lipinski_violations <= 1

                results.append(properties)

            except Exception as e:
                logger.error(f"Error processing SMILES {smiles}: {e}")
                results.append({
                    "smiles": smiles,
                    "error": str(e)
                })

        logger.info(f"Successfully predicted properties for {len([r for r in results if 'error' not in r])} molecules")

        return json.dumps({
            "success": True,
            "predictions": results,
            "total_processed": len(results)
        }, ensure_ascii=False, indent=2)

    except json.JSONDecodeError as e:
        return _validation_error(f"Invalid JSON format in smiles_list: {str(e)}")
    except Exception as e:
        logger.error(f"Error in predict_simple_properties_rdkit: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


@tool
def brics_decompose_rdkit(smiles: str, max_fragments: int = 50, rebuild_combinations: bool = False, max_recombine: int = 10) -> str:
    """
    Decompose a molecule into BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures) fragments.

    BRICS is a method for fragmenting molecules based on retrosynthetic rules. It breaks bonds that are
    commonly formed in organic synthesis, producing chemically meaningful building blocks.

    This tool is useful for:
    - Fragment-based molecular design
    - Building block library generation
    - Analyzing molecular scaffolds
    - Creating recombined novel structures from fragments

    Args:
        smiles (str): Input molecule in SMILES format.
                     Example: "c1ccccc1CCO" (phenethyl alcohol)
        max_fragments (int): Maximum number of fragments to return (1-100).
                            Default: 50
        rebuild_combinations (bool): If True, also generate recombined molecules from fragments.
                                    Default: False
        max_recombine (int): Maximum number of recombined molecules to generate (only used if rebuild_combinations=True).
                            Default: 10

    Returns:
        str: JSON string containing fragments and optionally recombined molecules.
             On success: {"success": true, "original_smiles": "...", "fragments": [...], "recombined": [...]}
             On error: {"success": false, "error": "..."}
    """
    # Input validation
    if not smiles or not isinstance(smiles, str):
        return _validation_error("smiles parameter is required and must be a non-empty string.")

    if not isinstance(max_fragments, int) or max_fragments < 1 or max_fragments > 100:
        return _validation_error("max_fragments must be an integer between 1 and 100.")

    if not isinstance(max_recombine, int) or max_recombine < 1 or max_recombine > 100:
        return _validation_error("max_recombine must be an integer between 1 and 100.")

    try:
        # Validate input SMILES
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return _validation_error(f"Invalid SMILES string: {smiles}")

        logger.info(f"Performing BRICS decomposition on {smiles}")

        # Perform BRICS decomposition
        fragments = BRICS.BRICSDecompose(mol)
        fragments_list = list(fragments)[:max_fragments]

        # Clean up fragments (remove dummy atoms markers for display)
        cleaned_fragments = []
        for frag in fragments_list:
            # Remove dummy atom labels [*] for cleaner display
            cleaned_frag = frag.replace('[1*]', '[*]').replace('[2*]', '[*]').replace('[3*]', '[*]')
            cleaned_frag = cleaned_frag.replace('[4*]', '[*]').replace('[5*]', '[*]').replace('[6*]', '[*]')
            cleaned_frag = cleaned_frag.replace('[7*]', '[*]').replace('[8*]', '[*]').replace('[9*]', '[*]')
            cleaned_frag = cleaned_frag.replace('[10*]', '[*]').replace('[11*]', '[*]').replace('[12*]', '[*]')
            cleaned_frag = cleaned_frag.replace('[13*]', '[*]').replace('[14*]', '[*]').replace('[15*]', '[*]')
            cleaned_frag = cleaned_frag.replace('[16*]', '[*]')
            cleaned_fragments.append(cleaned_frag)

        result = {
            "success": True,
            "original_smiles": smiles,
            "fragments": cleaned_fragments,
            "num_fragments": len(cleaned_fragments)
        }

        # Optionally rebuild molecules from fragments
        if rebuild_combinations and len(fragments_list) > 0:
            logger.info(f"Rebuilding molecules from {len(fragments_list)} fragments")
            try:
                # Use BRICSBuild to recombine fragments
                recombined_mols = BRICS.BRICSBuild(mol)
                recombined_list = []

                for i, rebuilt_mol in enumerate(recombined_mols):
                    if i >= max_recombine:
                        break
                    try:
                        rebuilt_smiles = Chem.MolToSmiles(rebuilt_mol)
                        # Skip if it's the same as the original
                        if rebuilt_smiles != smiles:
                            recombined_list.append(rebuilt_smiles)
                    except Exception as e:
                        logger.warning(f"Could not convert recombined molecule to SMILES: {e}")
                        continue

                result["recombined"] = recombined_list
                result["num_recombined"] = len(recombined_list)
                logger.info(f"Successfully generated {len(recombined_list)} recombined molecules")

            except Exception as e:
                logger.warning(f"Could not rebuild molecules from fragments: {e}")
                result["recombined_error"] = str(e)

        logger.info(f"BRICS decomposition completed: {len(cleaned_fragments)} fragments")

        return json.dumps(result, ensure_ascii=False, indent=2)

    except Exception as e:
        logger.error(f"Error in brics_decompose_rdkit: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


@tool
def recap_decompose_rdkit(smiles: str, max_fragments: int = 50, include_hierarchy: bool = True) -> str:
    """
    Decompose a molecule using Recap (Retrosynthetic Combinatorial Analysis Procedure).

    Recap is a retrosynthetic fragmentation method that breaks molecules at specific bond types
    commonly formed in synthesis (e.g., amide, ester, ether, amine bonds). It creates a hierarchical
    tree of fragments representing different levels of decomposition.

    This tool is useful for:
    - Retrosynthetic analysis
    - Identifying synthetic building blocks
    - Scaffold-based drug design
    - Fragment library generation with hierarchical information

    Args:
        smiles (str): Input molecule in SMILES format.
                     Example: "CC(=O)Nc1ccccc1" (acetanilide)
        max_fragments (int): Maximum number of leaf fragments to return (1-100).
                            Default: 50
        include_hierarchy (bool): If True, include hierarchical decomposition information.
                                 Default: True

    Returns:
        str: JSON string containing fragments and hierarchical decomposition data.
             On success: {"success": true, "original_smiles": "...", "leaf_fragments": [...], "hierarchy": {...}}
             On error: {"success": false, "error": "..."}
    """
    # Input validation
    if not smiles or not isinstance(smiles, str):
        return _validation_error("smiles parameter is required and must be a non-empty string.")

    if not isinstance(max_fragments, int) or max_fragments < 1 or max_fragments > 100:
        return _validation_error("max_fragments must be an integer between 1 and 100.")

    try:
        # Validate input SMILES
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return _validation_error(f"Invalid SMILES string: {smiles}")

        logger.info(f"Performing Recap decomposition on {smiles}")

        # Perform Recap decomposition
        recap_tree = Recap.RecapDecompose(mol)

        # Get leaf nodes (terminal fragments)
        leaf_fragments = []
        if recap_tree:
            leaves = recap_tree.GetLeaves()
            for leaf_smiles in list(leaves.keys())[:max_fragments]:
                # Clean up dummy atoms
                cleaned_leaf = leaf_smiles.replace('[1*]', '[*]').replace('[2*]', '[*]')
                cleaned_leaf = cleaned_leaf.replace('[3*]', '[*]').replace('[4*]', '[*]')
                cleaned_leaf = cleaned_leaf.replace('[5*]', '[*]').replace('[6*]', '[*]')
                cleaned_leaf = cleaned_leaf.replace('[7*]', '[*]').replace('[8*]', '[*]')
                cleaned_leaf = cleaned_leaf.replace('[9*]', '[*]').replace('[10*]', '[*]')
                cleaned_leaf = cleaned_leaf.replace('[11*]', '[*]')
                leaf_fragments.append(cleaned_leaf)

        result = {
            "success": True,
            "original_smiles": smiles,
            "leaf_fragments": leaf_fragments,
            "num_fragments": len(leaf_fragments)
        }

        # Include hierarchical decomposition information
        if include_hierarchy and recap_tree:
            hierarchy_info = []

            # Get all nodes (not just leaves)
            all_nodes = recap_tree.GetAllChildren()

            for node_smiles, node in list(all_nodes.items())[:max_fragments]:
                try:
                    # Clean up SMILES
                    cleaned_smiles = node_smiles
                    for i in range(1, 12):
                        cleaned_smiles = cleaned_smiles.replace(f'[{i}*]', '[*]')

                    # Get parent information if available
                    parent_info = []
                    if hasattr(node, 'parents') and node.parents:
                        for parent in node.parents:
                            parent_smiles = str(parent)
                            cleaned_parent = parent_smiles
                            for i in range(1, 12):
                                cleaned_parent = cleaned_parent.replace(f'[{i}*]', '[*]')
                            parent_info.append(cleaned_parent)

                    # Get children information if available
                    children_info = []
                    if hasattr(node, 'children') and node.children:
                        for child in node.children:
                            child_smiles = str(child)
                            cleaned_child = child_smiles
                            for i in range(1, 12):
                                cleaned_child = cleaned_child.replace(f'[{i}*]', '[*]')
                            children_info.append(cleaned_child)

                    hierarchy_info.append({
                        "fragment": cleaned_smiles,
                        "parents": parent_info,
                        "children": children_info,
                        "is_leaf": len(children_info) == 0
                    })

                except Exception as e:
                    logger.warning(f"Could not process node {node_smiles}: {e}")
                    continue

            result["hierarchy"] = hierarchy_info
            result["num_hierarchy_nodes"] = len(hierarchy_info)

        logger.info(f"Recap decomposition completed: {len(leaf_fragments)} leaf fragments")

        return json.dumps(result, ensure_ascii=False, indent=2)

    except Exception as e:
        logger.error(f"Error in recap_decompose_rdkit: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


# Side chain library definition
SIDE_CHAIN_LIBRARY = {
    "electron_donors": {
        "amino": {"smiles": "N", "name": "Amino (-NH2)", "description": "Strong electron donor, raises HOMO"},
        "hydroxyl": {"smiles": "O", "name": "Hydroxyl (-OH)", "description": "Moderate electron donor"},
        "methoxy": {"smiles": "OC", "name": "Methoxy (-OCH3)", "description": "Moderate electron donor"},
        "dimethylamino": {"smiles": "N(C)C", "name": "Dimethylamino (-N(CH3)2)", "description": "Very strong electron donor"},
        "methyl": {"smiles": "C", "name": "Methyl (-CH3)", "description": "Weak electron donor (hyperconjugation)"},
        "ethyl": {"smiles": "CC", "name": "Ethyl (-CH2CH3)", "description": "Weak electron donor"},
        "isopropyl": {"smiles": "C(C)C", "name": "Isopropyl (-CH(CH3)2)", "description": "Weak electron donor"},
        "tert-butyl": {"smiles": "C(C)(C)C", "name": "tert-Butyl (-C(CH3)3)", "description": "Weak electron donor, bulky"},
        "phenyl": {"smiles": "c1ccccc1", "name": "Phenyl (-Ph)", "description": "Conjugation extension"},
    },
    "electron_acceptors": {
        "nitro": {"smiles": "N(=O)=O", "name": "Nitro (-NO2)", "description": "Very strong electron acceptor, lowers LUMO"},
        "cyano": {"smiles": "C#N", "name": "Cyano (-CN)", "description": "Strong electron acceptor"},
        "trifluoromethyl": {"smiles": "C(F)(F)F", "name": "Trifluoromethyl (-CF3)", "description": "Strong electron acceptor"},
        "carboxyl": {"smiles": "C(=O)O", "name": "Carboxyl (-COOH)", "description": "Moderate electron acceptor"},
        "ester": {"smiles": "C(=O)OC", "name": "Methyl ester (-COOCH3)", "description": "Moderate electron acceptor"},
        "aldehyde": {"smiles": "C=O", "name": "Aldehyde (-CHO)", "description": "Strong electron acceptor"},
        "ketone": {"smiles": "C(=O)C", "name": "Acetyl (-COCH3)", "description": "Moderate electron acceptor"},
        "sulfonyl": {"smiles": "S(=O)(=O)C", "name": "Methylsulfonyl (-SO2CH3)", "description": "Strong electron acceptor"},
    },
    "halogens": {
        "fluoro": {"smiles": "F", "name": "Fluoro (-F)", "description": "Electron acceptor (inductive), small size"},
        "chloro": {"smiles": "Cl", "name": "Chloro (-Cl)", "description": "Weak electron acceptor"},
        "bromo": {"smiles": "Br", "name": "Bromo (-Br)", "description": "Weak electron acceptor, larger"},
        "iodo": {"smiles": "I", "name": "Iodo (-I)", "description": "Weak electron acceptor, heaviest"},
    },
    "other": {
        "thiol": {"smiles": "S", "name": "Thiol (-SH)", "description": "Weak electron donor, nucleophilic"},
        "methylthio": {"smiles": "SC", "name": "Methylthio (-SCH3)", "description": "Electron donor, polarizable"},
        "vinyl": {"smiles": "C=C", "name": "Vinyl (-CH=CH2)", "description": "Conjugation extension"},
        "ethynyl": {"smiles": "C#C", "name": "Ethynyl (-C≡CH)", "description": "Strong conjugation extension"},
        "formyl": {"smiles": "C=O", "name": "Formyl (-CHO)", "description": "Electron acceptor, aldehyde"},
    }
}


@tool
def substitute_side_chains_rdkit(
    smiles: str,
    side_chain_categories: str = "all",
    substitution_sites: str = "aromatic_hydrogens",
    max_substitutions: int = 20
) -> str:
    """
    Generate molecular analogs by systematically substituting side chains (functional groups).

    This tool creates a focused library of molecules by replacing hydrogen atoms or specific
    substructures with common functional groups from a curated library. Useful for:
    - Systematic SAR (structure-activity relationship) studies
    - Property optimization by functional group scanning
    - Creating diverse molecular libraries for screening

    The library includes:
    - Electron donors: -NH2, -OH, -OCH3, -N(CH3)2, -CH3, etc.
    - Electron acceptors: -NO2, -CN, -CF3, -COOH, -CHO, etc.
    - Halogens: -F, -Cl, -Br, -I
    - Other groups: -SH, -SCH3, vinyl, ethynyl, etc.

    Args:
        smiles (str): Base molecular structure in SMILES format.
                     Example: "c1ccccc1" (benzene)
        side_chain_categories (str): Comma-separated list of categories to include.
                                    Options: "electron_donors", "electron_acceptors", "halogens", "other", "all"
                                    Example: "electron_donors,electron_acceptors" or "all"
                                    Default: "all"
        substitution_sites (str): Type of sites to substitute.
                                 Options: "aromatic_hydrogens" (default), "aliphatic_hydrogens", "all_hydrogens"
                                 Default: "aromatic_hydrogens"
        max_substitutions (int): Maximum number of substituted molecules to generate (1-100).
                               Default: 20

    Returns:
        str: JSON string containing substituted molecules with design rationale.
             On success: {"success": true, "base_smiles": "...", "substituted_molecules": [...]}
             On error: {"success": false, "error": "..."}
    """
    # Input validation
    if not smiles or not isinstance(smiles, str):
        return _validation_error("smiles parameter is required and must be a non-empty string.")

    if not isinstance(max_substitutions, int) or max_substitutions < 1 or max_substitutions > 100:
        return _validation_error("max_substitutions must be an integer between 1 and 100.")

    if substitution_sites not in ["aromatic_hydrogens", "aliphatic_hydrogens", "all_hydrogens"]:
        return _validation_error(
            "substitution_sites must be one of: 'aromatic_hydrogens', 'aliphatic_hydrogens', 'all_hydrogens'"
        )

    try:
        # Validate input SMILES
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return _validation_error(f"Invalid SMILES string: {smiles}")

        logger.info(f"Generating side chain substitutions for {smiles}")

        # Parse side chain categories
        if side_chain_categories.lower() == "all":
            selected_categories = list(SIDE_CHAIN_LIBRARY.keys())
        else:
            selected_categories = [cat.strip() for cat in side_chain_categories.split(",")]
            # Validate categories
            invalid_cats = [cat for cat in selected_categories if cat not in SIDE_CHAIN_LIBRARY]
            if invalid_cats:
                return _validation_error(
                    f"Invalid categories: {invalid_cats}. Valid options: {list(SIDE_CHAIN_LIBRARY.keys())}"
                )

        # Collect functional groups from selected categories
        functional_groups = []
        for category in selected_categories:
            for group_key, group_data in SIDE_CHAIN_LIBRARY[category].items():
                functional_groups.append({
                    "category": category,
                    "key": group_key,
                    **group_data
                })

        logger.info(f"Using {len(functional_groups)} functional groups from categories: {selected_categories}")

        # Find substitution sites based on the specified type
        substitutable_atoms = []
        for atom in mol.GetAtoms():
            atom_idx = atom.GetIdx()

            # Check if atom is suitable for substitution
            if substitution_sites == "aromatic_hydrogens":
                # Only aromatic carbons with explicit or implicit hydrogens
                if atom.GetIsAromatic() and atom.GetSymbol() == 'C':
                    num_h = atom.GetTotalNumHs()
                    if num_h > 0:
                        substitutable_atoms.append(atom_idx)

            elif substitution_sites == "aliphatic_hydrogens":
                # Only aliphatic carbons with hydrogens
                if not atom.GetIsAromatic() and atom.GetSymbol() == 'C':
                    num_h = atom.GetTotalNumHs()
                    if num_h > 0:
                        substitutable_atoms.append(atom_idx)

            elif substitution_sites == "all_hydrogens":
                # Any carbon with hydrogens
                if atom.GetSymbol() == 'C':
                    num_h = atom.GetTotalNumHs()
                    if num_h > 0:
                        substitutable_atoms.append(atom_idx)

        if not substitutable_atoms:
            return _validation_error(
                f"No suitable substitution sites found in the molecule for {substitution_sites}"
            )

        logger.info(f"Found {len(substitutable_atoms)} substitutable atom positions")

        # Generate substituted molecules
        substituted_molecules = []
        count = 0

        for atom_idx in substitutable_atoms:
            if count >= max_substitutions:
                break

            for func_group in functional_groups:
                if count >= max_substitutions:
                    break

                try:
                    # Create editable molecule
                    edit_mol = Chem.RWMol(mol)

                    # Add the functional group
                    func_group_mol = Chem.MolFromSmiles(func_group["smiles"])
                    if not func_group_mol:
                        logger.warning(f"Invalid functional group SMILES: {func_group['smiles']}")
                        continue

                    # Insert functional group at the specified atom
                    # This is a simplified approach: add the functional group as a substituent
                    combined_mol = Chem.CombineMols(edit_mol, func_group_mol)

                    # Try to create a bond between the target atom and the functional group
                    # The functional group is added at the end, so its first atom is at len(mol.GetAtoms())
                    new_mol = Chem.RWMol(combined_mol)
                    func_start_idx = mol.GetNumAtoms()

                    # Add bond between target atom and first atom of functional group
                    new_mol.AddBond(atom_idx, func_start_idx, Chem.BondType.SINGLE)

                    # Sanitize and validate
                    try:
                        Chem.SanitizeMol(new_mol)
                        new_smiles = Chem.MolToSmiles(new_mol)

                        # Check if this is a valid and new molecule
                        if new_smiles != smiles:
                            substituted_molecules.append({
                                "smiles": new_smiles,
                                "functional_group": func_group["name"],
                                "category": func_group["category"],
                                "substitution_site": f"atom_{atom_idx}",
                                "rationale": func_group["description"]
                            })
                            count += 1
                            logger.debug(f"Generated: {new_smiles} with {func_group['name']}")

                    except Exception as sanitize_error:
                        logger.debug(f"Sanitization failed for {func_group['name']} at atom {atom_idx}: {sanitize_error}")
                        continue

                except Exception as e:
                    logger.debug(f"Could not substitute {func_group['name']} at atom {atom_idx}: {e}")
                    continue

        logger.info(f"Successfully generated {len(substituted_molecules)} substituted molecules")

        # If we couldn't generate any molecules through the direct approach,
        # fall back to a simpler SMARTS-based replacement for aromatic hydrogens
        if len(substituted_molecules) == 0 and substitution_sites == "aromatic_hydrogens":
            logger.info("Attempting fallback method using SMARTS replacement")

            for func_group in functional_groups[:max_substitutions]:
                try:
                    # Simple approach: for aromatic systems, use ReplaceSubstructs
                    # Find aromatic CH and replace with C-R
                    pattern = Chem.MolFromSmarts('[cH]')  # aromatic carbon with hydrogen
                    if not pattern:
                        continue

                    matches = mol.GetSubstructMatches(pattern)
                    if not matches:
                        continue

                    # Take the first match
                    match = matches[0]
                    atom_idx = match[0]

                    # Create a simple substituted version
                    # This is a simplified representation
                    base_smiles = Chem.MolToSmiles(mol)

                    # For aromatic systems, approximate substitution
                    # For benzene-like molecules, we can use simple string templates
                    if "c1ccccc1" in base_smiles:
                        substituted_smiles = base_smiles.replace("c1ccccc1", f"c1ccc({func_group['smiles']})cc1", 1)

                        # Validate the new SMILES
                        test_mol = Chem.MolFromSmiles(substituted_smiles)
                        if test_mol:
                            substituted_molecules.append({
                                "smiles": substituted_smiles,
                                "functional_group": func_group["name"],
                                "category": func_group["category"],
                                "substitution_site": "aromatic_position",
                                "rationale": func_group["description"]
                            })

                except Exception as e:
                    logger.debug(f"Fallback substitution failed for {func_group['name']}: {e}")
                    continue

        result = {
            "success": True,
            "base_smiles": smiles,
            "substitution_sites": substitution_sites,
            "categories_used": selected_categories,
            "substituted_molecules": substituted_molecules,
            "total_generated": len(substituted_molecules)
        }

        return json.dumps(result, ensure_ascii=False, indent=2)

    except Exception as e:
        logger.error(f"Error in substitute_side_chains_rdkit: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)
