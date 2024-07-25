from rdkit.Chem.rdmolfiles import MolFromPDBFile
from rdkit.Chem.rdchem import Mol

import numpy as np
import rdkit.Chem as Chem
from rdkit.Chem import AssignStereochemistry, HybridizationType, ChiralType, BondStereo, MolFromMol2File
import os
import sys
sys.path.append("../../")
from src.utils.constants import ES_THRESHOLD, METAL_OX_STATES, ES_DIST_THRESHOLD, VDW_RADII, PDB_DATA
from src.utils.debugging import time_block

def get_vdw_radius(symbol):
    return VDW_RADII.get(symbol, 200) / 100  # 200 is Default value for unknown elements, convert to Angstrom

def pdb_to_rdkit_mol(pdb_filepath: str, sanitize: bool = True): 

    #check if the file exists
    if not os.path.exists(pdb_filepath):
        raise FileNotFoundError(f"File {pdb_filepath} not found")

    mol = MolFromPDBFile(pdb_filepath, removeHs=True, sanitize = sanitize)

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

            #check if the atom is a hydrogen 
            if parts[5] == 'H':
                continue            
            
            # Extract charge from the line (assuming it's the last column)
            if len(parts) >= 9:  # Ensure there are enough columns
                try:
                    charge = float(parts[-1])
                    charges.append(charge)
                except ValueError:
                    # Handle cases where charge cannot be converted to float
                    charges.append(None)  # Or handle differently based on your needs

    return charges

def get_covalent_bond_edge_features(mol) -> tuple:

    #loop over all the bonds in the molecule and get the edge features
    edge_features, edge_indices = [], []

    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()

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
        bond_features = (bond_order, aromaticity, conjugation, ring, stereochemistry)

        edge_indices.append((atom1, atom2))
        edge_indices.append((atom2, atom1))

        edge_features.append(bond_features)
        edge_features.append(bond_features)

    return np.array(edge_features, dtype='float64'), np.array(edge_indices, dtype='int64')

def vdw_interactions(mol, conf, atom_1, atom_2):
    """
    Determine if there is a van der Waals interaction between atom 1 and atom 2 based on whether the interatomic distance is less 
    than ES_DIST_THRESHOLD.
    
    Parameters:
    mol (Chem.Mol): The RDKit molecule.
    atom_1 (int): The index of the first atom.
    atom_2 (int): The index of the second atom.
    
    Returns:
    tuple: (0, 0, 0, 0, 0) if there is a van der Waals interaction, None otherwise 
    """

    pos_1 = conf.GetAtomPosition(atom_1)
    pos_2 = conf.GetAtomPosition(atom_2)
    distance = pos_1.Distance(pos_2)

    # Get the van der Waals radii of the two atoms
    atom_1_symbol = mol.GetAtomWithIdx(atom_1).GetSymbol()
    atom_2_symbol = mol.GetAtomWithIdx(atom_2).GetSymbol()
    vdw_radius_1 = get_vdw_radius(atom_1_symbol)
    vdw_radius_2 = get_vdw_radius(atom_2_symbol)

    # Check if the distance is less than or equal to the sum of the van der Waals radii
    if distance <= (vdw_radius_1 + vdw_radius_2):
        return (0, 0, 0, 0, 0)

    return None

def get_vdw_radius(symbol: str) -> float:
    # Placeholder function to get van der Waals radius for a given atom symbol
    vdw_radii = {
        'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.80, 'S': 1.80, 'Cl': 1.75,
        # Add more elements as needed
    }
    return vdw_radii.get(symbol, 1.5)  # Default value if the element is not in the dictionary

def get_vdw_interaction_edge_features(mol) -> tuple:
    """
    Extracts edge features for van der Waals interactions from a given RDKit molecule object.

    Parameters:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object.

    Returns:
        edge_features (np.ndarray): A 2D array of shape (num_edges, num_edge_features) containing edge features.
        edge_indices (np.ndarray): A 2D array of shape (num_edges, 2) containing the indices of the atoms connected by each edge.
    """
    # Get the conformer
    conf = mol.GetConformer()

    # Get positions of all atoms
    positions = conf.GetPositions()
    
    # Compute pairwise distances using matrix operations
    dist_matrix = np.sqrt(np.sum((positions[:, np.newaxis, :] - positions[np.newaxis, :, :]) ** 2, axis=-1))

    # Get van der Waals radii for all atoms
    vdw_radii = np.array([get_vdw_radius(mol.GetAtomWithIdx(i).GetSymbol()) for i in range(mol.GetNumAtoms())])

    # Create a matrix of summed van der Waals radii
    vdw_sum_matrix = vdw_radii[:, np.newaxis] + vdw_radii[np.newaxis, :]

    # Find pairs of atoms within the van der Waals interaction distance
    atom_pairs = np.transpose(np.where((dist_matrix <= vdw_sum_matrix) & (dist_matrix > 0)))

    # Initialize lists to store edge features and indices
    edge_indices, edge_features = [], []

    # Iterate over atom pairs
    for pair in atom_pairs:
        atom1, atom2 = pair
        
        # Append edge indices to list, duplicating to account for both directions
        edge_indices.append((atom1, atom2))

        # Append edge features to list, duplicating to account for both directions
        edge_feature = (0, 0, 0, 0, 0)
        edge_features.append(edge_feature)
    
    return np.array(edge_features, dtype='float64'), np.array(edge_indices, dtype='int64')

def ionic_interaction_using_pcs(mol: Chem.Mol, conf, atom_1: int, atom_2: int):
    """
    Determine if there is an electrostatic interaction between atom 1 and atom 2 based on whether the coulombic interaction
    between the two atoms is less than or equal to some threshold imported from src/utils/constants.py.

    Parameters:
    mol (Chem.Mol): The RDKit molecule. NOTE: We assume that the molecule has _TriposAtomCharges charges assigned to the atoms.
    atom_1 (int): The index of the first atom.
    atom_2 (int): The index of the second atom.
    
    Returns:
    tuple: (0, 0, 0, 0, 0) if there is an electrostatic interaction, None otherwise 

    Conversion Factor Calculation:
    The Coulombic interaction energy between two point charges q1 and q2 separated by a distance r in vacuum is given by:

        E = k * q1 * q2 / r

    Where:
    - E is the energy in Joules (J).
    - k is Coulomb's constant, 8.987 x 10^9 N m^2 C^-2.
    - q1 and q2 are the charges in Coulombs (C).
    - r is the distance in meters (m).

    To convert this energy into kJ/mol, the following conversions are applied:
    - 1 e (elementary charge) = 1.602 x 10^-19 C
    - Avogadro's number (N_A) = 6.022 x 10^23 mol^-1
    - 1 Angstrom (Å) = 10^-10 meters (m)

    Therefore, the conversion factor from (charge_1 * charge_2) / distance to kJ/mol is:

        Conversion factor = (k * N_A * (1.602 x 10^-19)^2) * 10^10

    Plugging in the values:

        Conversion factor = (8.987 x 10^9 * 6.022 x 10^23 * (1.602 x 10^-19)^2) * 10^10
                          ≈ 1388.93 kJ mol^-1 Å^-1
    """
    # Constants
    CONVERSION_FACTOR_KJ = 1388.93  # kJ mol^-1 Å^-1

    #Ensure that the atoms are not covanlently bonded
    if mol.GetBondBetweenAtoms(atom_1, atom_2):
        return None

    # Get atomic charges using Gasteiger charges or another suitable method
    charge_1 = float(mol.GetAtomWithIdx(atom_1).GetProp('_TriposAtomCharges'))
    charge_2 = float(mol.GetAtomWithIdx(atom_2).GetProp('_TriposAtomCharges'))
    
    pos_1 = conf.GetAtomPosition(atom_1)
    pos_2 = conf.GetAtomPosition(atom_2)
    distance = pos_1.Distance(pos_2)

    coulombic_interaction_kj = (charge_1 * charge_2 / distance) * CONVERSION_FACTOR_KJ

    # Check if the Coulombic interaction is less than or equal to the threshold
    if abs(coulombic_interaction_kj) >= ES_THRESHOLD:
        return (0, 0, 0, 0, 0, coulombic_interaction_kj)
    else:
        return None
    
def ionic_interaction_using_distance(mol: Chem.Mol, conf, atom_1: int, atom_2: int):
    """
    Determine if there is an electrostatic interaction between atom 1 and atom 2 based on whether the distance between the two atoms is less than 
    ES_DIST_THRESHOLD (defined in src/utils/constants.py).
    
    Parameters:
    mol (Chem.Mol): The RDKit molecule. NOTE: We assume that the molecule has _TriposAtomCharges charges assigned to the atoms.
    atom_1 (int): The index of the first atom.
    atom_2 (int): The index of the second atom.
    
    Returns:
    tuple: (0, 0, 0, 0, distance) if there is an electrostatic interaction, None otherwise 

    """ 
    
    pos_1 = conf.GetAtomPosition(atom_1)
    pos_2 = conf.GetAtomPosition(atom_2)
    distance = pos_1.Distance(pos_2)

    # Check if the Coulombic interaction is less than or equal to the threshold
    if distance <= ES_DIST_THRESHOLD:
        return (0, 0, 0, 0, 0, distance)
    
    else:

        return None

def get_ionic_interaction_edge_features(mol) -> tuple:
    """
    Extracts edge features for ionic interactions from a given RDKit molecule object.

    Parameters:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object.

    Returns:
        edge_features (np.ndarray): A 2D array of shape (num_edges, num_edge_features) containing edge features.
        edge_indices (np.ndarray): A 2D array of shape (num_edges, 2) containing the indices of the atoms connected by each edge.
    """
    # Get the conformer
    conf = mol.GetConformer()

    # Get positions of all atoms
    positions = conf.GetPositions()
    
    # Compute pairwise distances using matrix operations
    dist_matrix = np.sqrt(np.sum((positions[:, np.newaxis, :] - positions[np.newaxis, :, :]) ** 2, axis=-1))
    
    # Find pairs of atoms within the distance threshold
    atom_pairs = np.transpose(np.where((dist_matrix <= ES_DIST_THRESHOLD) & (dist_matrix > 0)))
    
    # Initialize lists to store edge features and indices
    edge_indices, edge_features = [], []

    # Iterate over atom pairs
    for pair in atom_pairs:
        atom1, atom2 = pair
        if atom1 < atom2: 
            pos_1 = positions[atom1]
            pos_2 = positions[atom2]
            distance = np.linalg.norm(pos_1 - pos_2)
            
            # Append edge indices to list, duplicating to account for both directions
            edge_indices.append((atom1, atom2))
            edge_indices.append((atom2, atom1))

            # Append edge features to list, duplicating to account for both directions
            edge_feature = (0, 0, 0, 0, distance)
            edge_features.append(edge_feature)
            edge_features.append(edge_feature)
    
    return np.array(edge_features, dtype='float64'), np.array(edge_indices, dtype='int64')

def add_charges_to_molecule(mol, charges):

    if len(mol.GetAtoms()) != len(charges):
        raise ValueError(f"The number of charges {len(charges)} does not match the number of atoms {len(mol.GetAtoms())} in the molecule.")
    
    for atom, charge in zip(mol.GetAtoms(), charges):
        #check if the atom is a metal 
        if atom.GetSymbol() in METAL_OX_STATES.keys():
            print(atom.GetSymbol())
            charge = METAL_OX_STATES[atom.GetSymbol()] 
        
        atom.SetDoubleProp('_TriposAtomCharges', charge)

def get_node_features(mol: Mol) -> np.ndarray:
    """
    Extracts node features from a given RDKit molecule object.
    Parameters:

        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object. It should have a "protein_or_ligand" property for each atom

    Returns:
        node_features (np.ndarray): A 2D array of shape (num_atoms, num_node_features) containing node features:
            (protein_or_ligand_id, atomic_num, atomic_mass, aromatic_indicator, ring_indicator, hybridization,
        chirality, num_heteroatoms, degree, num_hydrogens, formal_charge, num_radical_electrons)
    """
    
    node_features  = []

    # Iterate over each atom in the molecule and calculate node features
    for atom in mol.GetAtoms():

        protein_or_ligand_id = [-1 if atom.GetProp('protein_or_ligand') == 'protein' else 1][0]
        
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
        num_hydrogens = len([bond for bond in atom.GetBonds() if bond.GetOtherAtom(atom).GetAtomicNum() == 1])

        formal_charge = atom.GetFormalCharge()
        num_radical_electrons = atom.GetNumRadicalElectrons()

        # Append node features to list
        node_features.append((protein_or_ligand_id, atomic_num, atomic_mass, aromatic_indicator, ring_indicator, hybridization,
        chirality, num_heteroatoms, degree, num_hydrogens, formal_charge, num_radical_electrons))

    return np.array(node_features, dtype='float64')

def rdkit_mol_featurizer(mol: Mol,
                         pairwise_function: str,
                         ) -> tuple: 
    """
    Featurizes an RDKit molecule by extracting node features, edge features, and edge indices.
    
    Args:
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule to be featurized.
        pairwise_function (str): The type of pairwise function to be used for generating edge features.
            Supported values are "ionic interactions", "covalent bonds", and "vdw interactions".
    
    Returns:
        tuple: A tuple containing the node features, edge features, and edge indices as NumPy arrays.
    """
    try: 
        # for atom in mol.GetAtoms():
        #     atom.UpdatePropertyCache()
        AssignStereochemistry(mol, force=True, cleanIt=True)
    except Exception as e:
        #add to the error log file 
        with open(os.path.join(PDB_DATA, 'featurization_error.log'), 'a') as file:
            file.write(f"Error processing {mol.GetProp('_Name')}: {e}")

    node_features = get_node_features(mol)

    if pairwise_function == "ionic interactions":
        edge_features, edge_indices = get_ionic_interaction_edge_features(mol)
    elif pairwise_function == "covalent bonds":
        edge_features, edge_indices = get_covalent_bond_edge_features(mol)
    elif pairwise_function == "vdw interactions":
        edge_features, edge_indices = get_vdw_interaction_edge_features(mol)
    else:
        raise ValueError(f"Pairwise function {pairwise_function} not supported")

    return np.array(node_features, dtype='float64'), np.array(edge_features, dtype='float64'), np.array(edge_indices, dtype='int64')