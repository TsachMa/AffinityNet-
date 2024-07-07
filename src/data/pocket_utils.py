from openbabel import openbabel
from openbabel import pybel
import numpy as np
from rdkit.Chem import Mol

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

def find_pocket_atoms_RDKit(mol : Mol,
                            protein_or_ligand_ids: list,
                            num_atoms_in_protein):

    #identify the pocket atoms
    #first, add all of the ligand atoms 
    pocket_atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if protein_or_ligand_ids[atom.GetIdx()] == 1]
   
    #add the protein atoms that are within 4 angstroms of the ligand atoms
    for l_atom in list(mol.GetAtoms())[num_atoms_in_protein:]: #loop over all of the ligand atoms 
        for p_atom in list(mol.GetAtoms())[:num_atoms_in_protein]: #loop over all of the protein atoms 
            distance = np.linalg.norm(get_atom_coordinates(mol)[l_atom.GetIdx()] - get_atom_coordinates(mol)[p_atom.GetIdx()])

            if distance < 4.0: # 4 angstroms 
                #if the index is not already in the pocket atom indices, add it
                if p_atom.GetIdx() not in pocket_atom_indices:
                    pocket_atom_indices.append(p_atom.GetIdx())
    return pocket_atom_indices

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
