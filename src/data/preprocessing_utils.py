from openbabel import openbabel, pybel
import os
from tqdm import tqdm 
import sys
import re

def remove_water_molecules(input_file, use_CLI = True):
    """
    Args: 
        input_file (str): The path to the input file (mol2)
        output_file (str): The path to the output file (mol2)
    """
    water_residue_names = ["HOH", "H2O", "TIP3"]  # Add more if needed
    output_file = input_file.replace(".mol2", "_no_water.mol2")

    if use_CLI:
        command = f"obabel {input_file} -O {output_file}"
        for residue_name in water_residue_names:
            command += f" -xr {residue_name}"
        
        log_file = "error.log"
        try:
            os.system(command)
        except Exception as e:
            with open(log_file, 'a') as f:
                f.write(f"Error processing file: {input_file}\n")

    else:
        try: 
            # Read the molecule
            mol = next(pybel.readfile("mol2", input_file))
        except StopIteration:
            print(f"Error reading file {input_file}")
            sys.exit(1)

        # Identify water residues
        waters_to_remove = []
        for residue in mol.residues:
            for water_residue_name in water_residue_names:
                if water_residue_name in residue.OBResidue.GetName():
                    waters_to_remove.append(residue)

        # Remove water residues
        for water in waters_to_remove:
            mol.OBMol.DeleteResidue(water.OBResidue)

        # Write the molecule to the output file
        mol.write("mol2", output_file, overwrite=True)

    

def combine_protein_ligand(protein_file, ligand_file, output_file):
    """
    Combines a protein and ligand mol2 file into a single protein-ligand complex mol2 file.

    Args:
        protein_file (str): The path to the protein mol2 file.
        ligand_file (str): The path to the ligand mol2 file.
        output_file (str): The path to the output mol2 file.
    """
    # Read the protein and ligand molecules
    protein = next(pybel.readfile("mol2", protein_file))
    ligand = next(pybel.readfile("mol2", ligand_file))

    # Combine the molecules
    combined_mol = pybel.Molecule(protein.OBMol)
    combined_mol.OBMol += ligand.OBMol

    # Write the combined molecule to the output file
    combined_mol.write("mol2", output_file, overwrite=True)

def add_charges_and_save(input_filepath):
    """
    Args:
        input_filepath (str): The path to the input file (mol2)
    """

    # Extract pdb id
    base_name = input_filepath.split("/")[-1].split("_")[0]
    output_file = f"{base_name}_complex_charged.mol2"
    # Initialize OBConversion
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "mol2")

    # Read the protein and ligand molecules
    complex = pybel.readfile("mol2", input_filepath).__next__().OBMol
    
    # Calculate charges
    chargeModel = openbabel.OBChargeModel.FindType("mmff94")
    if not chargeModel.ComputeCharges(complex):
        raise ValueError("Failed to compute charges for the molecule ", input_filepath)
    
    # Write the output file
    obConversion.WriteFile(complex, output_file)


def get_filenames(directory, regex, recursive=True):
    """
    Retrieves filenames from a directory matching a specified regular expression.

    Args:
        directory (str): The path to the directory.
        regex (str): The regular expression pattern to match filenames.
        recursive (bool): If True, search directories recursively. Default is True.

    Returns:
        List[str]: A list of filenames that match the regex pattern.
    """
    matched_files = []
    pattern = re.compile(regex)

    if recursive:
        for root, _, files in os.walk(directory):
            for file in files:
                if pattern.match(file):
                    matched_files.append(os.path.join(root, file))
    else:
        for file in os.listdir(directory):
            if os.path.isfile(os.path.join(directory, file)) and pattern.match(file):
                matched_files.append(os.path.join(directory, file))

    return matched_files