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

from services.pubchem_service import PubChemService
from services.exceptions import NotFoundError, ValidationError, ServiceError

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
        base_mol = Chem.MolFromSmiles(base_smiles)
        if not base_mol:
            return _validation_error(f"Invalid SMILES string: {base_smiles}")

        logger.info(f"Generating analogs for {base_smiles} with strategy: {design_strategy}")

        # Get canonical SMILES for base molecule
        canonical_base = Chem.MolToSmiles(base_mol)

        # Collect all candidates
        candidates = []

        # Functional groups to add based on strategy
        conjugation_groups = [
            ("C=C", "vinyl", "Added vinyl group to extend conjugation"),
            ("C#C", "ethynyl", "Added ethynyl group for extended conjugation"),
            ("c1ccccc1", "phenyl", "Added phenyl ring to extend conjugation"),
        ]

        push_pull_pairs = [
            ("N", "C#N", "amino (donor)", "cyano (acceptor)", "Push-pull: Amino (donor) and cyano (acceptor)"),
            ("N(C)C", "C#N", "dimethylamino (donor)", "cyano (acceptor)", "Push-pull: Dimethylamino (strong donor) and cyano (acceptor)"),
            ("O", "N(=O)=O", "hydroxyl (donor)", "nitro (acceptor)", "Push-pull: Hydroxyl (donor) and nitro (acceptor)"),
            ("N", "C(F)(F)F", "amino (donor)", "trifluoromethyl (acceptor)", "Push-pull: Amino (donor) and trifluoromethyl (acceptor)"),
            ("OC", "C(=O)O", "methoxy (donor)", "carboxyl (acceptor)", "Push-pull: Methoxy (donor) and carboxyl (acceptor)"),
        ]

        donor_groups = [
            ("N", "amino", "Added amino group (strong electron donor, raises HOMO)"),
            ("O", "hydroxyl", "Added hydroxyl group (moderate electron donor)"),
            ("OC", "methoxy", "Added methoxy group (moderate electron donor)"),
            ("C", "methyl", "Added methyl group (weak electron donor via hyperconjugation)"),
        ]

        acceptor_groups = [
            ("C#N", "cyano", "Added cyano group (strong electron acceptor, lowers LUMO)"),
            ("N(=O)=O", "nitro", "Added nitro group (very strong electron acceptor)"),
            ("C(F)(F)F", "trifluoromethyl", "Added trifluoromethyl group (strong electron acceptor)"),
            ("C(=O)C", "acetyl", "Added acetyl group (moderate electron acceptor)"),
            ("C(=O)O", "carboxyl", "Added carboxyl group (moderate electron acceptor)"),
        ]

        # Strategy-based generation
        if design_strategy == "conjugation" or design_strategy == "diverse":
            # Add conjugation-extending groups
            for group_smiles, group_name, rationale in conjugation_groups:
                analog = _try_add_group_to_aromatic_h(base_mol, group_smiles, rationale)
                if analog:
                    candidates.append(analog)
                    if len(candidates) >= num_to_generate:
                        break

        if design_strategy == "push_pull" or design_strategy == "diverse":
            # Add push-pull substituents
            for donor, acceptor, donor_name, acceptor_name, rationale in push_pull_pairs:
                analog = _try_add_push_pull_substituents(base_mol, donor, acceptor, rationale)
                if analog:
                    candidates.append(analog)
                    if len(candidates) >= num_to_generate:
                        break

        if design_strategy == "heteroatom" or design_strategy == "diverse":
            # Replace aromatic CH with heteroatoms
            heteroatom_analogs = _try_heteroatom_substitution(base_mol, num_to_generate)
            candidates.extend(heteroatom_analogs)
            if len(candidates) >= num_to_generate:
                candidates = candidates[:num_to_generate]

        if design_strategy == "diverse":
            # Add electron donors
            for group_smiles, group_name, rationale in donor_groups:
                if len(candidates) >= num_to_generate:
                    break
                analog = _try_add_group_to_aromatic_h(base_mol, group_smiles, rationale)
                if analog:
                    candidates.append(analog)

            # Add electron acceptors
            for group_smiles, group_name, rationale in acceptor_groups:
                if len(candidates) >= num_to_generate:
                    break
                analog = _try_add_group_to_aromatic_h(base_mol, group_smiles, rationale)
                if analog:
                    candidates.append(analog)

        # Deduplicate by SMILES
        seen_smiles = {canonical_base}  # Exclude base molecule
        validated_candidates = []
        for candidate in candidates:
            if candidate["smiles"] not in seen_smiles:
                seen_smiles.add(candidate["smiles"])
                validated_candidates.append(candidate)
                if len(validated_candidates) >= num_to_generate:
                    break

        logger.info(f"Successfully generated {len(validated_candidates)} valid analogs")

        return json.dumps({
            "success": True,
            "base_smiles": canonical_base,
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


def _try_add_group_to_aromatic_h(base_mol: Chem.Mol, group_smiles: str, rationale: str) -> Optional[Dict[str, str]]:
    """
    Try to add a functional group to an aromatic hydrogen position.

    Args:
        base_mol: Base molecule
        group_smiles: SMILES of the group to add
        rationale: Design rationale

    Returns:
        Dict with 'smiles' and 'rationale' if successful, None otherwise
    """
    try:
        # Find aromatic carbons with hydrogens
        aromatic_ch_pattern = Chem.MolFromSmarts('[cH]')
        if not aromatic_ch_pattern:
            return None

        matches = base_mol.GetSubstructMatches(aromatic_ch_pattern)
        if not matches:
            # Try aliphatic carbons if no aromatic CH found
            aliphatic_ch_pattern = Chem.MolFromSmarts('[CH]')
            if aliphatic_ch_pattern:
                matches = base_mol.GetSubstructMatches(aliphatic_ch_pattern)

        if not matches:
            return None

        # Try to substitute at the first match
        atom_idx = matches[0][0]

        # Create editable copy
        edit_mol = Chem.RWMol(base_mol)

        # Parse functional group
        group_mol = Chem.MolFromSmiles(group_smiles)
        if not group_mol:
            return None

        # Add functional group atoms to molecule
        group_start_idx = edit_mol.GetNumAtoms()
        for atom in group_mol.GetAtoms():
            edit_mol.AddAtom(Chem.Atom(atom.GetAtomicNum()))

        # Add bonds within the functional group
        for bond in group_mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx() + group_start_idx
            end_idx = bond.GetEndAtomIdx() + group_start_idx
            edit_mol.AddBond(begin_idx, end_idx, bond.GetBondType())

        # Connect functional group to base molecule
        edit_mol.AddBond(atom_idx, group_start_idx, Chem.BondType.SINGLE)

        # Sanitize and validate
        try:
            Chem.SanitizeMol(edit_mol)
            new_smiles = Chem.MolToSmiles(edit_mol)
            return {"smiles": new_smiles, "rationale": rationale}
        except Exception:
            return None

    except Exception as e:
        logger.debug(f"Failed to add group {group_smiles}: {e}")
        return None


def _try_add_push_pull_substituents(base_mol: Chem.Mol, donor_smiles: str, acceptor_smiles: str, rationale: str) -> Optional[Dict[str, str]]:
    """
    Try to add push-pull substituents (donor and acceptor) to aromatic positions.

    Args:
        base_mol: Base molecule
        donor_smiles: SMILES of electron donor group
        acceptor_smiles: SMILES of electron acceptor group
        rationale: Design rationale

    Returns:
        Dict with 'smiles' and 'rationale' if successful, None otherwise
    """
    try:
        # Find aromatic carbons with hydrogens
        aromatic_ch_pattern = Chem.MolFromSmarts('[cH]')
        if not aromatic_ch_pattern:
            return None

        matches = base_mol.GetSubstructMatches(aromatic_ch_pattern)
        if len(matches) < 2:
            return None

        # Try to add at two different positions (ideally para or meta)
        atom_idx_1 = matches[0][0]
        atom_idx_2 = matches[min(len(matches) - 1, 2)][0]  # Try to get some separation

        # Create editable copy
        edit_mol = Chem.RWMol(base_mol)

        # Add donor group
        donor_mol = Chem.MolFromSmiles(donor_smiles)
        if not donor_mol:
            return None

        donor_start_idx = edit_mol.GetNumAtoms()
        for atom in donor_mol.GetAtoms():
            edit_mol.AddAtom(Chem.Atom(atom.GetAtomicNum()))

        for bond in donor_mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx() + donor_start_idx
            end_idx = bond.GetEndAtomIdx() + donor_start_idx
            edit_mol.AddBond(begin_idx, end_idx, bond.GetBondType())

        edit_mol.AddBond(atom_idx_1, donor_start_idx, Chem.BondType.SINGLE)

        # Add acceptor group
        acceptor_mol = Chem.MolFromSmiles(acceptor_smiles)
        if not acceptor_mol:
            return None

        acceptor_start_idx = edit_mol.GetNumAtoms()
        for atom in acceptor_mol.GetAtoms():
            edit_mol.AddAtom(Chem.Atom(atom.GetAtomicNum()))

        for bond in acceptor_mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx() + acceptor_start_idx
            end_idx = bond.GetEndAtomIdx() + acceptor_start_idx
            edit_mol.AddBond(begin_idx, end_idx, bond.GetBondType())

        edit_mol.AddBond(atom_idx_2, acceptor_start_idx, Chem.BondType.SINGLE)

        # Sanitize and validate
        try:
            Chem.SanitizeMol(edit_mol)
            new_smiles = Chem.MolToSmiles(edit_mol)
            return {"smiles": new_smiles, "rationale": rationale}
        except Exception:
            return None

    except Exception as e:
        logger.debug(f"Failed to add push-pull groups: {e}")
        return None


def _try_heteroatom_substitution(base_mol: Chem.Mol, max_candidates: int = 5) -> List[Dict[str, str]]:
    """
    Try to substitute aromatic CH with heteroatoms (N, O, S).

    Args:
        base_mol: Base molecule
        max_candidates: Maximum number of candidates to generate

    Returns:
        List of dicts with 'smiles' and 'rationale'
    """
    candidates = []

    try:
        # Find aromatic carbons
        aromatic_c_pattern = Chem.MolFromSmarts('[c]')
        if not aromatic_c_pattern:
            return candidates

        matches = base_mol.GetSubstructMatches(aromatic_c_pattern)
        if not matches:
            return candidates

        # Try substituting with different heteroatoms
        heteroatoms = [
            (7, "N", "nitrogen"),  # Nitrogen
            (8, "O", "oxygen"),    # Oxygen
            (16, "S", "sulfur"),   # Sulfur
        ]

        for atom_num, symbol, name in heteroatoms:
            if len(candidates) >= max_candidates:
                break

            # Try substituting the first aromatic carbon
            atom_idx = matches[0][0]

            # Create editable copy
            edit_mol = Chem.RWMol(base_mol)

            # Get the carbon atom
            carbon_atom = edit_mol.GetAtomWithIdx(atom_idx)

            # Check if substitution makes sense (e.g., not already heteroatom)
            if carbon_atom.GetAtomicNum() == 6:  # Is carbon
                # Replace with heteroatom
                carbon_atom.SetAtomicNum(atom_num)
                carbon_atom.SetFormalCharge(0)

                # Try to sanitize
                try:
                    Chem.SanitizeMol(edit_mol)
                    new_smiles = Chem.MolToSmiles(edit_mol)
                    candidates.append({
                        "smiles": new_smiles,
                        "rationale": f"Replaced aromatic carbon with {name} to create heteroaromatic system"
                    })
                except Exception:
                    logger.debug(f"Failed to substitute with {name}")
                    continue

    except Exception as e:
        logger.debug(f"Error in heteroatom substitution: {e}")

    return candidates


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

        # Get canonical base SMILES
        canonical_base = Chem.MolToSmiles(mol)

        # Generate substituted molecules using the helper function
        substituted_molecules = []
        seen_smiles = {canonical_base}  # Exclude base molecule

        for func_group in functional_groups:
            if len(substituted_molecules) >= max_substitutions:
                break

            # Build rationale
            rationale = f"{func_group['description']} (category: {func_group['category']})"

            # Try to add the functional group
            result = _try_add_group_to_aromatic_h(mol, func_group["smiles"], rationale)

            if result and result["smiles"] not in seen_smiles:
                # Add additional metadata
                substituted_molecules.append({
                    "smiles": result["smiles"],
                    "functional_group": func_group["name"],
                    "category": func_group["category"],
                    "rationale": result["rationale"]
                })
                seen_smiles.add(result["smiles"])
                logger.debug(f"Generated: {result['smiles']} with {func_group['name']}")

        logger.info(f"Successfully generated {len(substituted_molecules)} substituted molecules")

        result = {
            "success": True,
            "base_smiles": canonical_base,
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


@tool
def search_compound_smiles_pubchem(compound_name_en: str) -> str:
    """
    Search PubChem for a compound by its English name and retrieve the canonical SMILES representation.

    **IMPORTANT**: This tool ONLY accepts compound names in English. If the user provides a
    compound name in another language (e.g., Japanese, Spanish), you MUST translate it to
    English before calling this tool.

    This tool is critical for ensuring accurate molecular structure representation when
    users specify compound names or molecular scaffolds. Always use this tool when:
    - A user requests a molecule by name (e.g., "azulene", "naphthalene", "benzene")
    - You need to verify the correct SMILES for a molecular scaffold
    - You are unsure about the SMILES representation of a compound

    **Examples:**
    - User says "アズレン" → Translate to "azulene" → Call this tool with "azulene"
    - User says "ナフタレン" → Translate to "naphthalene" → Call this tool with "naphthalene"
    - User says "benzene" → Directly call this tool with "benzene"

    Args:
        compound_name_en (str): The compound name in English.
                               Example: "azulene", "naphthalene", "benzene"
                               MUST be in English.

    Returns:
        str: JSON string containing the canonical SMILES and compound information.
             On success: {
                 "success": true,
                 "smiles": "...",
                 "compound_info": {
                     "cid": 123,
                     "iupac_name": "...",
                     "molecular_formula": "...",
                     "molecular_weight": 123.45
                 }
             }
             On error: {"success": false, "error": "..."}
    """
    # Input validation
    if not compound_name_en or not isinstance(compound_name_en, str):
        return _validation_error("compound_name_en parameter is required and must be a non-empty string.")

    if not compound_name_en.strip():
        return _validation_error("compound_name_en cannot be empty or whitespace only.")

    try:
        logger.info(f"Searching PubChem for compound: {compound_name_en}")

        # Use PubChemService to retrieve SMILES
        pubchem_service = PubChemService()
        result = pubchem_service.get_compound_smiles(compound_name_en, search_type='name')

        logger.info(f"Successfully retrieved SMILES for '{compound_name_en}': {result['smiles']}")

        return json.dumps({
            "success": True,
            "smiles": result['smiles'],
            "compound_info": result['compound_info']
        }, ensure_ascii=False, indent=2)

    except NotFoundError as e:
        logger.warning(f"Compound not found: {compound_name_en}")
        return json.dumps({
            "success": False,
            "error": f"Compound not found in PubChem: {compound_name_en}. Please check the spelling or try a different name."
        }, ensure_ascii=False, indent=2)
    except ValidationError as e:
        logger.error(f"Validation error for compound search: {e}")
        return _validation_error(str(e))
    except ServiceError as e:
        logger.error(f"PubChem service error: {e}")
        return json.dumps({
            "success": False,
            "error": f"PubChem service error: {str(e)}"
        }, ensure_ascii=False, indent=2)
    except Exception as e:
        logger.error(f"Error in search_compound_smiles_pubchem: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)
