from openbabel import openbabel
from openbabel import pybel
import numpy as np
from rdkit.Chem import Mol
from tqdm import tqdm 

from rdkit import Chem
from rdkit.Chem import AllChem

import sys
sys.path.append('../../')
from src.utils.constants import POCKET_THRESHOLD
from src.utils.debugging import time_block

def read_molecule(filename):
    """Read a molecule from a file."""
    return next(pybel.readfile(filename.split('.')[-1], filename))

def find_pocket_atoms_OB(protein, 
                      ligand, 
                      radius):
    """Find protein atoms within a certain radius of the ligand
     using OpenBabel."""
    pocket_atoms = []
    for protein_atom in protein:
        for ligand_atom in ligand:
            distance = protein_atom.OBAtom.GetDistance(ligand_atom.OBAtom)
            if distance <= radius:
                pocket_atoms.append(protein_atom)
                break
    return pocket_atoms

def get_atom_coordinates(mol):
    """
    Extracts atom coordinates from a molecule.
    """
    conformer = mol.GetConformer()
    coords = []
    for atom in mol.GetAtoms():
        pos = conformer.GetAtomPosition(atom.GetIdx())
        coords.append((pos.x, pos.y, pos.z))
    return np.array(coords)

def remove_waters(mol):
    water_residue_names = {"HOH", "WAT"}
    atoms_to_remove = []

    for atom in mol.GetAtoms():
        residue_info = atom.GetPDBResidueInfo()
        if residue_info and residue_info.GetResidueName().strip() in water_residue_names:
            atoms_to_remove.append(atom.GetIdx())
    
    editable_mol = Chem.EditableMol(mol)
    for atom_idx in reversed(atoms_to_remove):
        editable_mol.RemoveAtom(atom_idx)
    
    return editable_mol.GetMol()

def get_centroid(mol):
    conf = mol.GetConformer()
    positions = conf.GetPositions()
    centroid = np.mean(positions, axis=0)
    return centroid

def get_max_distance_from_centroid(mol, centroid):
    conf = mol.GetConformer()
    max_distance = 0.0
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        distance = np.linalg.norm(pos - centroid)
        if distance > max_distance:
            max_distance = distance
    return max_distance

def get_atoms_within_distance(mol, centroid, distance, by_residue=False):
    """
    Returns a list of atom indices within a specified distance from a given centroid.

    Args:
        mol (Chem.Mol): The molecule object.
        centroid (numpy.ndarray): The centroid position.
        distance (float): The maximum distance from the centroid.
        by_residue (bool, optional): If True, returns atom indices grouped by residue. Defaults to False.

    Returns:
        list: A list of atom indices within the specified distance.

    """
    conf = mol.GetConformer()

    atom_positions = np.array(conf.GetPositions())

    distances = np.linalg.norm(atom_positions - centroid, axis=1)
    within_distance_mask = distances <= distance

    atom_indices_within_distance = np.where(within_distance_mask)[0]

    if by_residue:
        residue_numbers = np.array([atom.GetPDBResidueInfo().GetResidueNumber() for atom in mol.GetAtoms()])

        chain_ids = np.array([atom.GetPDBResidueInfo().GetChainId() for atom in mol.GetAtoms()])
        
        # Get residues and chains within the distance
        residues_within_distance = set(zip(residue_numbers[within_distance_mask], chain_ids[within_distance_mask]))
        
        # Find atom indices that match the residue number and chain id
        atom_indices_within_distance = np.where(
            np.array([(residue_numbers[i], chain_ids[i]) in residues_within_distance for i in range(len(residue_numbers))])
        )[0]
        
    return list(atom_indices_within_distance)

def tag_atoms(mol, tag):
    for atom in mol.GetAtoms():
        atom.SetProp('protein_or_ligand', tag)

def combine_and_filter(ligand_mol, protein_mol, constant_distance = POCKET_THRESHOLD, by_residue=False):
    ligand_centroid = get_centroid(ligand_mol)

    max_ligand_distance = get_max_distance_from_centroid(ligand_mol, ligand_centroid)
    threshold_distance = max_ligand_distance + constant_distance

    atoms_to_include = get_atoms_within_distance(protein_mol, ligand_centroid, threshold_distance, by_residue)

    # Assuming 'protein_mol' is your RDKit molecule object
    atoms_to_include_set = set([int(atom_idx) for atom_idx in atoms_to_include])

    # Create a new editable molecule
    new_mol = Chem.EditableMol(Chem.Mol())

    # Map original atom indices to new atom indices
    atom_map = {}

    # Add the atoms to the new molecule along with their properties
    for atom_idx in atoms_to_include:
        atom = protein_mol.GetAtomWithIdx(int(atom_idx))
        new_idx = new_mol.AddAtom(atom)
        atom_map[atom_idx] = new_idx

    # Create a conformer to hold the coordinates
    conformer = Chem.Conformer(len(atoms_to_include_set))

    # Set the coordinates for each atom in the new molecule
    for atom_idx in atoms_to_include:
        pos = protein_mol.GetConformer().GetAtomPosition(int(atom_idx))
        conformer.SetAtomPosition(atom_map[atom_idx], pos)

    # Add the bonds to the new molecule
    for bond in protein_mol.GetBonds():
        start_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()
        if start_atom in atoms_to_include_set and end_atom in atoms_to_include_set:
            new_mol.AddBond(atom_map[start_atom], atom_map[end_atom], bond.GetBondType())

    # Get the final molecule with coordinates
    editable_protein_mol = new_mol.GetMol()

    #add the conformer
    editable_protein_mol.AddConformer(conformer)

    # Copy atom properties
    for atom_idx in atoms_to_include:
        atom = protein_mol.GetAtomWithIdx(int(atom_idx))
        new_atom = editable_protein_mol.GetAtomWithIdx(atom_map[atom_idx])
        for prop in atom.GetPropNames():
            new_atom.SetProp(prop, atom.GetProp(prop))

    # Copy molecule properties
    for prop in protein_mol.GetPropNames():
        editable_protein_mol.SetProp(prop, protein_mol.GetProp(prop))

    filtered_protein_mol = editable_protein_mol

    tag_atoms(ligand_mol, "ligand")
    tag_atoms(filtered_protein_mol, "protein")

    complex_mol = Chem.CombineMols(filtered_protein_mol, ligand_mol)

    return complex_mol

def identify_atom_sources(complex_mol):
    atom_sources = []
    for atom in complex_mol.GetAtoms():
        source = atom.GetProp('source')
        atom_sources.append(source)
    return atom_sources

def save_pocket_atoms_OB(pocket_atoms, output_filename):
    """Save the pocket atoms to a PDB file."""
    mol = pybel.Molecule(openbabel.OBMol())
    for atom in pocket_atoms:
        mol.OBMol.AddAtom(atom.OBAtom)
    mol.write("pdb", output_filename, overwrite=True)

def find_and_save_pocket_atoms_OB(protein_filename, ligand_filename, radius, output_filename):
    """Find and save protein atoms within a certain radius of the ligand.
    
    Args:   
        protein_filename (str): The filename of the protein PDB file.
        ligand_filename (str): The filename of the ligand PDB file.
        radius (float): The radius in angstroms.
        output_filename (str): The filename of the output PDB file
        
    Returns:
        None
    """

    protein = read_molecule(protein_filename)
    ligand = read_molecule(ligand_filename)
    pocket_atoms = find_pocket_atoms_OB(protein, ligand, radius)
    save_pocket_atoms_OB(pocket_atoms, output_filename)
    print(f"Pocket atoms saved to {output_filename}")
