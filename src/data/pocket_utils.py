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

def find_pocket_atoms_RDKit(mol: Chem.Mol, protein_or_ligand_ids: list, num_atoms_in_protein):
    # Identify ligand atoms
    ligand_atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if protein_or_ligand_ids[atom.GetIdx()] == 1]
    
    # Get coordinates of ligand atoms
    ligand_coords = get_atom_coordinates(mol)[num_atoms_in_protein:]
    
    # Compute the centroid of the ligand atoms
    ligand_centroid = np.mean(ligand_coords, axis=0)
    
    # Compute the maximum distance from the centroid to any ligand atom
    max_distance = max(np.linalg.norm(ligand_coords - ligand_centroid, axis=1)) + POCKET_THRESHOLD  # Adding 4 Å to the max distance

    pocket_atom_indices = set(ligand_atom_indices)

    # Get coordinates of all atoms
    all_coords = get_atom_coordinates(mol)
    
    # Loop over protein atoms and check if they are within the bounding sphere
    for p_atom in tqdm(mol.GetAtoms()):
        p_idx = p_atom.GetIdx()
        if p_idx >= num_atoms_in_protein:  # Skip ligand atoms
            continue
        
        distance_to_centroid = np.linalg.norm(all_coords[p_idx] - ligand_centroid)
        if distance_to_centroid <= max_distance:
            pocket_atom_indices.add(p_idx)
                    
    return list(pocket_atom_indices)

def save_pocket_atoms(pocket_atoms, output_filename):
    """Save the pocket atoms to a PDB file."""
    mol = pybel.Molecule(openbabel.OBMol())
    for atom in pocket_atoms:
        mol.OBMol.AddAtom(atom.OBAtom)
    mol.write("pdb", output_filename, overwrite=True)

def find_and_save_pocket_atoms(protein_filename, ligand_filename, radius, output_filename):
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
    save_pocket_atoms(pocket_atoms, output_filename)
    print(f"Pocket atoms saved to {output_filename}")
