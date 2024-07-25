from multiprocessing import Pool, cpu_count
import os
from tqdm import tqdm
from rdkit import Chem
import h5py
import sys 
sys.path.append('../')
from src.data.utils import mol2_to_rdkit_mol, pdb_to_rdkit_mol, rdkit_mol_featurizer
from src.data.pocket_utils import remove_waters, combine_and_filter
from src.utils.constants import PDB_DATA
import warnings

PARALLEL = False 
CLEAN = True 

def process_protein(protein_name):
    #check if the protein has already been processed
    if os.path.exists(os.path.join(PDB_DATA, protein_name, 'protein_data.h5')):
        return None

    try: 
        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always")

        ligand_mol2_path = os.path.join(PDB_DATA, f'{protein_name}/{protein_name}_ligand.mol2')
        protein_pdb_path = os.path.join(PDB_DATA, f'{protein_name}/{protein_name}_protein.pdb')

        ligand_mol = mol2_to_rdkit_mol(ligand_mol2_path, sanitize=False)
        protein_mol = pdb_to_rdkit_mol(protein_pdb_path, sanitize=False)

        pocket_mol_res = combine_and_filter(ligand_mol, remove_waters(protein_mol), by_residue=True)

        output_path = os.path.join(PDB_DATA, f'{protein_name}/{protein_name}_pocket_res.pdb')
        pdb_block = Chem.MolToPDBBlock(pocket_mol_res)
        
        with open(output_path, 'w') as file:
            file.write(pdb_block)

        vdw_graph = rdkit_mol_featurizer(pocket_mol_res, "vdw interactions")
        ionic_graph = rdkit_mol_featurizer(pocket_mol_res, "ionic interactions")
        cov_graph = rdkit_mol_featurizer(pocket_mol_res, "covalent bonds")

        protein_data = {}
        protein_data = extend_protein_data_dict(protein_data, protein_name, vdw_graph, ionic_graph, cov_graph)

        #cache out the protein data
        save_protein_data_to_hdf5(os.path.join(PDB_DATA, protein_name, 'protein_data.h5'), protein_data)

        if caught_warnings:
            with open(os.path.join(PDB_DATA, 'featurization_error.log'), 'a') as file:
                for warning in caught_warnings:
                    file.write(f'{protein_name}: Warning: {str(warning.message)}\n')

        return None
    
    except Exception as e:

        with open(os.path.join(PDB_DATA, 'featurization_error.log'), 'a') as file:
            file.write(f'{protein_name}: {str(e)}\n')

        return None

def extend_protein_data_dict(protein_data, protein_name, vdw_graph, ionic_graph, cov_graph):
    protein_data[protein_name] = {}
    protein_data[protein_name]['vdw'] = {}
    protein_data[protein_name]['vdw']['node_features'] = vdw_graph[0]
    protein_data[protein_name]['vdw']['edge_features'] = vdw_graph[1]
    protein_data[protein_name]['vdw']['edge_list'] = vdw_graph[2]

    protein_data[protein_name]['ionic'] = {}
    protein_data[protein_name]['ionic']['node_features'] = ionic_graph[0]
    protein_data[protein_name]['ionic']['edge_features'] = ionic_graph[1]
    protein_data[protein_name]['ionic']['edge_list'] = ionic_graph[2]

    protein_data[protein_name]['cov'] = {}
    protein_data[protein_name]['cov']['node_features'] = cov_graph[0]
    protein_data[protein_name]['cov']['edge_features'] = cov_graph[1]
    protein_data[protein_name]['cov']['edge_list'] = cov_graph[2]

    return protein_data

# Function to save data to HDF5
def save_protein_data_to_hdf5(file_name, protein_data):
    with h5py.File(file_name, 'w') as f:
        for protein, graphs in protein_data.items():
            protein_group = f.create_group(protein)
            for graph_name, graph_data in graphs.items():
                graph_group = protein_group.create_group(graph_name)
                for data_name, data in graph_data.items():
                    graph_group.create_dataset(data_name, data=data)

if __name__ == "__main__":
    pdb_data = PDB_DATA
    protein_data = {}

    if CLEAN:
        #clear the error log at PDB/featurization_error.log
        with open(os.path.join(PDB_DATA, 'featurization_error.log'), 'w') as file:
            file.write('')

    protein_dirs = sorted([d for d in os.listdir(pdb_data) if os.path.isdir(os.path.join(pdb_data, d))])

    if CLEAN:
        #loop through the protein directories and remove the protein_data.h5 file
        for protein in protein_dirs:
            if os.path.exists(os.path.join(PDB_DATA, protein, 'protein_data.h5')):
                os.remove(os.path.join(PDB_DATA, protein, 'protein_data.h5'))

    if PARALLEL:
        with Pool(cpu_count() // 2) as pool:
            results = list(tqdm(pool.imap(process_protein, protein_dirs), total=len(protein_dirs)))
    else:
        for protein in tqdm(protein_dirs):
            process_protein(protein)

    # for result in results:
    #     protein_name, vdw_graph, ionic_graph, cov_graph = result
    #     protein_data = extend_protein_data_dict(protein_data, protein_name, vdw_graph, ionic_graph, cov_graph)

    # save_protein_data_to_hdf5(os.path.join(pdb_data, 'protein_data.h5'), protein_data)
