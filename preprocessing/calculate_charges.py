import requests
import sys
import os
import argparse
from tqdm import tqdm 

def add_mol2_charges(complex_mol2_file):
    # Upload the pocket mol2 file to the ACC2 API
    r = requests.post('https://acc2-api.biodata.ceitec.cz/send_files', files={'file[]': open(complex_mol2_file, 'rb')})
    
    # Check if the upload was successful
    if r.status_code != 200:
        raise Exception(f"Failed to upload the file. Status code: {r.status_code}, Response: {r.text}")
    
    # Obtain ID number for uploaded file
    r_id = list(r.json()['structure_ids'].values())[0]
    
    # Calculate charges using eqeq method
    r_out = requests.get(f'https://acc2-api.biodata.ceitec.cz/calculate_charges?structure_id={r_id}&method=eqeq&generate_mol2=true')
    
    # Check if the charge calculation was successful
    if r_out.status_code != 200:
        raise Exception(f"Failed to calculate charges. Status code: {r_out.status_code}, Response: {r_out.text}")
    
    output_filename = complex_mol2_file.replace('.mol2', '_charged.mol2')
    # Save output mol2 file
    with open(output_filename, 'wb') as f:
        f.write(r_out.content)
    print(f"Charges added and saved to '{output_filename}'.")

def main():
    # #take as an argument the path to your data
    # parser = argparse.ArgumentParser(description='Create cleaned dataset')
    # parser.add_argument('--PDB_data_path', type=str, help='path to PDBbind dataset')

    # args = parser.parse_args()

    # PDBbind_data_path = os.path.join(args.PDB_data_path, 'PDBbind_2020_data.csv')
    # general_set_PDBs_path = os.path.join(args.PDB_data_path, 'v2020-other-PL')
    # directories = [PDBbind_data_path, general_set_PDBs_path]

    directories = ['/Users/tsachmackey/dfs/affinity_net/PDBbind/v2020-other-PL']

    for directory in directories:
        subdirs = [sub_dir for sub_dir in os.listdir(directory) if os.path.isdir(os.path.join(directory, sub_dir))]
        for sub_dir in tqdm(subdirs):
            for file in os.listdir(os.path.join(directory, sub_dir)):
                if file.endswith('complex.mol2'):
                    complex_mol2_file = os.path.join(directory, sub_dir, file)
                    add_mol2_charges(complex_mol2_file)

if __name__ == "__main__":
    main()