from rdkit.Chem.rdmolfiles import MolFromPDBFile
from rdkit.Chem.rdchem import Mol
import numpy as np
from rdkit.Chem import AddHs, AssignStereochemistry, HybridizationType, ChiralType, BondStereo, MolFromMol2File
from rdkit.Chem.AllChem import ComputeGasteigerCharges
import os
from pocket_utils import get_atom_coordinates, find_pocket_atoms_RDKit

def pdb_to_rdkit_mol(pdb_filepath: str): 

    #check if the file exists
    if not os.path.exists(pdb_filepath):
        raise FileNotFoundError(f"File {pdb_filepath} not found")

    mol = MolFromPDBFile(pdb_filepath, removeHs=True)

    return mol

def mol2_to_rdkit_mol(mol2_filepath: str, sanitize: bool = True): 

    #check if the file exists
    if not os.path.exists(mol2_filepath):
        raise FileNotFoundError(f"File {mol2_filepath} not found")

    mol = MolFromMol2File(mol2_filepath, sanitize=sanitize, removeHs=False)

    return mol

def get_protein_or_ligand_ids(pdb_filepath: str) -> list: 
    """
    Extracts whether an atom is a protein or a ligand from a given PDB file.
    Parameters:
        pdb_filepath (str): Path to the PDB file.
    Returns:
        atom_types (list): A list of strings containing the type of each atom in the PDB file.
        it is either -1 for protein or 1 for ligand.

    Description: 
        All atoms in lines starting with "ATOM" are considered as protein atoms and all atoms 
        in lines starting with "HETATM" are considered as ligand atoms.
    """

    atom_types = []
    with open(pdb_filepath, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_types.append(-1)
            elif line.startswith('HETATM'):
                atom_types.append(1)

    return atom_types

def extract_charges_from_mol2(file_path):
    charges = []

    # Open the MOL2 file for reading
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Flag to track if we are in the @<TRIPOS>ATOM section
    in_atom_section = False

    # Iterate through lines in the file
    for line in lines:
        line = line.strip()
        
        # Check for the start of @<TRIPOS>ATOM section
        if line.startswith('@<TRIPOS>ATOM'):
            in_atom_section = True
            continue
        
        # Check for the end of @<TRIPOS>ATOM section
        if in_atom_section and line.startswith('@<TRIPOS>'):
            break
        
        # Process lines within @<TRIPOS>ATOM section
        if in_atom_section:
            # Split line by whitespace
            parts = line.split()
            
            # Extract charge from the line (assuming it's the last column)
            if len(parts) >= 9:  # Ensure there are enough columns
                try:
                    charge = float(parts[-1])
                    charges.append(charge)
                except ValueError:
                    # Handle cases where charge cannot be converted to float
                    charges.append(None)  # Or handle differently based on your needs

    return charges

def get_node_features(mol: Mol,
                        protein_or_ligand_ids: list = None) -> np.ndarray:
    """
    Extracts node features from a given RDKit molecule object.
    Parameters:

        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object.
        protein_or_ligand_ids (list): A list of ints describing whether each atom in the molecule is a protein (-1) or a ligand (1) atom.
    Returns:
        node_features (np.ndarray): A 2D array of shape (num_atoms, num_node_features) containing node features:
            (protein_or_ligand_id, atomic_num, atomic_mass, aromatic_indicator, ring_indicator, hybridization, chirality, num_heteroatoms, degree, num_hydrogens, partial_charge, formal_charge, num_radical_electrons)
    """
    node_features  = []
    # Iterate over each atom in the molecule and calculate node features
    for atom in mol.GetAtoms():
        protein_or_ligand_id = protein_or_ligand_ids[atom.GetIdx()] 
        
        # Calculate node features
        atomic_num = atom.GetAtomicNum()
        atomic_mass = atom.GetMass()
        aromatic_indicator = int(atom.GetIsAromatic())
        ring_indicator = int(atom.IsInRing())
        hybridization_tag = atom.GetHybridization()
        if hybridization_tag == HybridizationType.SP:
            hybridization = 1
        elif hybridization_tag == HybridizationType.SP2:
            hybridization = 2
        elif hybridization_tag == HybridizationType.SP3:
            hybridization = 3
        else:
            hybridization = 0

        chiral_tag = atom.GetChiralTag()

        if chiral_tag == ChiralType.CHI_TETRAHEDRAL_CW:
            chirality = 1

        elif chiral_tag == ChiralType.CHI_TETRAHEDRAL_CCW:
            chirality = -1
        else:
            chirality = 0
        num_heteroatoms = len([bond for bond in atom.GetBonds() if bond.GetOtherAtom(atom).GetAtomicNum() != atom.GetAtomicNum()])
        degree = atom.GetDegree()
        num_hydrogens = atom.GetTotalNumHs()
        partial_charge = atom.GetDoubleProp('_TriposAtomCharges')
        formal_charge = atom.GetFormalCharge()
        num_radical_electrons = atom.GetNumRadicalElectrons()
        # Append node features to list
        node_features.append((protein_or_ligand_id, atomic_num, atomic_mass, aromatic_indicator, ring_indicator, hybridization,
        chirality, num_heteroatoms, degree, num_hydrogens, partial_charge, formal_charge, num_radical_electrons))

    return np.array(node_features, dtype='float64')

def covalent_bonds(mol, atom_1, atom_2):
    """
    Determine if there is a covalent bond between two specified atoms in a molecule and get bond features.
    
    Parameters:
    mol (Chem.Mol): The RDKit molecule.
    atom_1 (int): The index of the first atom.
    atom_2 (int): The index of the second atom.
    
    Returns:
    tuple: A tuple containing bond order, aromaticity, conjugation, ring, and stereochemistry.
           If there is no bond, returns (0, 0, 0, 0, 0).
    """
    bond = mol.GetBondBetweenAtoms(atom_1, atom_2)

    if bond: 
        # Calculate edge features
        bond_order = bond.GetBondTypeAsDouble()
        aromaticity = int(bond.GetIsAromatic())
        conjugation = int(bond.GetIsConjugated())
        ring = int(bond.IsInRing())
        stereochemistry_tag = bond.GetStereo()
        
        if stereochemistry_tag == BondStereo.STEREOZ:
            stereochemistry = 1
        elif stereochemistry_tag == BondStereo.STEREOE:
            stereochemistry = -1
        else:
            stereochemistry = 0
    else:
        bond_order, aromaticity, conjugation, ring, stereochemistry = 0, 0, 0, 0, 0

    return (bond_order, aromaticity, conjugation, ring, stereochemistry)


def get_edge_features(mol: Mol,
                      pocket_atom_indices: list,
                      pairwise_function: callable,
                      ) -> tuple:
    """

    Extracts edge features from a given RDKit molecule object.
    Parameters:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object.
        pocket_atom_indices (list): A list of ints containing the indices of the atoms in the pocket.
        pairwise_function (callable): A function that takes two atoms and returns edge features.


    Returns:
        edge_features (np.ndarray): A 2D array of shape (num_bonds, num_edge_features) containing edge features.
        edge_indices (np.ndarray): A 2D array of shape (num_bonds, 2) containing the indices of the atoms connected by each bond.
    
    """

    # Initialize a list to store edge features and indices
    edge_indices, edge_features = [], []

    #for every atom in the pocket, create an edge between atoms if they are within 4 angstroms of each other
    for i, atom1 in enumerate(pocket_atom_indices):
        atom_i = mol.GetAtomWithIdx(atom1)

        for j, atom2 in enumerate(pocket_atom_indices):
            atom_j = mol.GetAtomWithIdx(atom2)

            if j > i: #only consider the upper triangle of the matrix
                i_j_edge_features = pairwise_function(mol, atom1, atom2)

                if i_j_edge_features:
                    # Append edge indices to list, duplicating to account for both directions
                    edge_indices.append((atom1, atom2))
                    edge_indices.append((atom2, atom1))

                    # Append edge features to list, duplicating to account for both directions
                    edge_features.append(i_j_edge_features)
                    edge_features.append(i_j_edge_features)

    return np.array(edge_features, dtype='float64'), np.array(edge_indices, dtype='int64')

def add_charges_to_molecule(mol, charges):
    if len(mol.GetAtoms()) != len(charges):
        raise ValueError(f"The number of charges {len(charges)} does not match the number of atoms {len(mol.GetAtoms())} in the molecule.")
    
    for atom, charge in zip(mol.GetAtoms(), charges):
        atom.SetDoubleProp('_TriposAtomCharges', charge)

def rdkit_mol_featurizer(mol: Mol, 
                         charges: list, 
                         protein_or_ligand_ids: list = None, 
                         num_atoms_in_protein: int = None,
                         non_convalent_edges: bool = True,
                         ) -> tuple: 
    """
    Extracts graph features from a given RDKit molecule object.
    Parameters:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object.
        protein_or_ligand_ids (list): A list of ints describing whether each atom in the molecule is a protein (-1) or a ligand (1) atom. 
    Returns:
        node_features (np.ndarray): A 2D array of shape (num_atoms, num_node_features) containing node features:
            (protein_or_ligand_id, atomic_num, atomic_mass, aromatic_indicator, ring_indicator, hybridization, chirality, num_heteroatoms, degree, num_hydrogens, partial_charge, formal_charge, num_radical_electrons)
        edge_features (np.ndarray): A 2D array of shape (num_bonds, num_edge_features) containing edge features.
        edge_indices (np.ndarray): A 2D array of shape (num_bonds, 2) containing the indices of the atoms connected by each bond.
    """

    AssignStereochemistry(mol, force=True, cleanIt=True)

    if not protein_or_ligand_ids:
        protein_or_ligand_ids = [-1 if atom.GetIdx() < num_atoms_in_protein else 1 for atom in mol.GetAtoms()]

    add_charges_to_molecule(mol, charges)

    node_features = get_node_features(mol, protein_or_ligand_ids)
    
    pocket_atom_indices = find_pocket_atoms_RDKit(mol, protein_or_ligand_ids, num_atoms_in_protein)
    pocket_node_features = [node_features[i] for i in pocket_atom_indices]

    if non_convalent_edges:
        edge_features, edge_indices = get_edge_features(mol, pocket_atom_indices)
    else:
        edge_features, edge_indices = None, None

    return np.array(pocket_node_features, dtype='float64'), np.array(edge_features, dtype='float64'), np.array(edge_indices, dtype='int64')