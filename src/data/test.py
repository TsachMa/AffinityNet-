import numpy as np
import unittest
import sys
sys.path.append("../../")

from src.data.utils import get_node_features, rdkit_mol_featurizer, mol2_to_rdkit_mol, get_protein_or_ligand_ids, extract_charges_from_mol2
from src.data.pocket_utils import read_molecule, find_pocket_atoms_OB, find_and_save_pocket_atoms
from src.utils.constants import PROJECT_ROOT
HOME_DIR = PROJECT_ROOT

class TestRDKitMolFeaturizer(unittest.TestCase):
    def setUp(self):
        self.mol = mol2_to_rdkit_mol(f"{HOME_DIR}/test_data/pdb/test_methionine.mol2")
        self.protein_or_ligand_id = [-1, -1, -1, -1, -1, -1, -1, -1]
        self.charges = extract_charges_from_mol2(f"{HOME_DIR}/test_data/pdb/test_methionine.mol2")
        self.node_features, self.edge_features, self.edge_indices = rdkit_mol_featurizer(self.mol, self.charges, self.protein_or_ligand_id, non_convalent_edges=False)
    
    def test_charge_extractor(self):
        expected_charges = [-1.04, -0.074, 0.329, -0.545, -0.397, -0.458, 0.141, -0.661, 0.323, 0.321, 0.334, 0.162, 0.215, 0.192, 0.175, 0.197, 0.199, 0.194, 0.195, 0.198]
        np.testing.assert_array_almost_equal(self.charges, expected_charges)

    def test_rdkit_mol_featurizer(self):
        #we know that this receptor is just the methionine molecule
        #the features are: 
        # (atomic_num, atomic_mass, aromatic_indicator, ring_indicator, hybridization, chirality, 
        # num_heteroatoms, degree, num_hydrogens, partial_charge, formal_charge, num_radical_electrons)
        expected_node_feature_for_met_N = (-1, #protein_or_ligand_id
                                           7, #atomic_num
                                 14.007, #atomic_mass
                                 0,  #aromatic_indicator
                                 0,  #ring_indicator
                                 3,  #hybridization
                                 0,  #chirality
                                 1,  #num_heteroatoms
                                 1,  #degree
                                 2,  #num_hydrogens
                                 -3.21724048e-01,  #partial_charge
                                 0,  #formal_charge
                                 0) #num_radical_electrons
        
        expected_node_feature_for_alpha_C = (-1, #protein_or_ligand_id
                                             6, #atomic_num
                                    12.011, #atomic_mass
                                    0,  #aromatic_indicator
                                    0,  #ring_indicator
                                    3,  #hybridization
                                    -1,  #chirality
                                    1,  #num_heteroatoms
                                    3,  #degree
                                    1,  #num_hydrogens
                                    6.12758477e-02,  #partial_charge
                                    0,  #formal_charge
                                    0)
        
        expected_node_feature_for_met_ketone_C = (-1, #protein_or_ligand_id
                                                  6, #atomic_num
                                    12.011, #atomic_mass
                                    0,  #aromatic_indicator
                                    0,  #ring_indicator
                                    2,  #hybridization
                                    0,  #chirality
                                    1,  #num_heteroatoms
                                    2,  #degree
                                    1,  #num_hydrogens
                                    1.36206953e-01,  #partial_charge
                                    0,  #formal_charge
                                    0)
        
        expected_node_feature_for_met_ketone_O = (-1, #protein_or_ligand_id
                                                  8, #atomic_num
                                    15.999, #atomic_mass
                                    0,  #aromatic_indicator
                                    0,  #ring_indicator
                                    2,  #hybridization
                                    0,  #chirality
                                    1,  #num_heteroatoms
                                    1,  #degree
                                    0,  #num_hydrogens
                                    -3.01642718e-01,  #partial_charge
                                    0,  #formal_charge
                                    0)
        
        expected_node_feature_for_met_sidechain_C1 = (-1, #protein_or_ligand_id
                                                      6, #atomic_num
                                    12.011, #atomic_mass
                                    0,  #aromatic_indicator
                                    0,  #ring_indicator
                                    3,  #hybridization
                                    0,  #chirality
                                    0,  #num_heteroatoms
                                    2,  #degree
                                    2,  #num_hydrogens
                                    -2.20069145e-02,  #partial_charge
                                    0,  #formal_charge
                                    0)
        
        expected_node_feature_for_met_sidechain_C2 = (-1,
                                                      6, #atomic_num
                                    12.011, #atomic_mass
                                    0,  #aromatic_indicator
                                    0,  #ring_indicator
                                    3,  #hybridization
                                    0,  #chirality
                                    1,  #num_heteroatoms
                                    2,  #degree
                                    2,  #num_hydrogens
                                    -5.12891136e-03,  #partial_charge
                                    0,  #formal_charge
                                    0)
        
        expected_node_feature_for_met_sidechain_S = (-1,  #protein_or_ligand_id
                                                     16, #atomic_num
                                    32.067, #atomic_mass
                                    0,  #aromatic_indicator
                                    0,  #ring_indicator
                                    3,  #hybridization
                                    0,  #chirality
                                    2,  #num_heteroatoms
                                    2,  #degree
                                    0,  #num_hydrogens
                                    -1.65359663e-01,  #partial_charge
                                    0,  #formal_charge
                                    0)
        
        expected_node_feature_for_met_sidechain_C3 = (-1, #protein_or_ligand_id
                                                      6, #atomic_num
                                    12.011, #atomic_mass
                                    0,  #aromatic_indicator
                                    0,  #ring_indicator
                                    3,  #hybridization
                                    0,  #chirality
                                    1,  #num_heteroatoms
                                    1,  #degree
                                    3,  #num_hydrogens
                                    -1.83989336e-02,  #partial_charge
                                    0,  #formal_charge
                                    0)
        
        #the bond features are: 
        # (bond_order, aromaticity, conjugation, ring, stereochemistry)
        expected_bond_features = np.array([[1., 0., 0., 0., 0.],
                                           [1., 0., 0., 0., 0.],
                                           [1., 0., 0., 0., 0.],
                                           [1., 0., 0., 0., 0.],
                                           [2., 0., 0., 0., 0.],
                                           [2., 0., 0., 0., 0.],
                                           [1., 0., 0., 0., 0.],
                                           [1., 0., 0., 0., 0.],
                                           [1., 0., 0., 0., 0.],
                                           [1., 0., 0., 0., 0.],
                                           [1., 0., 0., 0., 0.],
                                           [1., 0., 0., 0., 0.],
                                           [1., 0., 0., 0., 0.],
                                           [1., 0., 0., 0., 0.]])
        
        #the edge indices are: 
        # (idx1, idx2)
        expected_edge_indices = np.array([[1, 0],
                                                    [0, 1],
                                                    [2, 1],
                                                    [1, 2],
                                                    [3, 2],
                                                    [2, 3],
                                                    [4, 1],
                                                    [1, 4],
                                                    [5, 4],
                                                    [4, 5],
                                                    [6, 5],
                                                    [5, 6],
                                                    [7, 6],
                                                    [6, 7]])
        
        #assuming that the ordering of the nodes is 
        # N, alpha_C, met_ketone_C, met_ketone_O, met_sidechain_C1, met_sidechain_C2, met_sidechain_S, met_sidechain_C3
        ordered_list_of_expected_node_features = [expected_node_feature_for_met_N, 
                                                  expected_node_feature_for_alpha_C,
                                                    expected_node_feature_for_met_ketone_C,
                                                    expected_node_feature_for_met_ketone_O,
                                                    expected_node_feature_for_met_sidechain_C1,
                                                    expected_node_feature_for_met_sidechain_C2,
                                                    expected_node_feature_for_met_sidechain_S,
                                                    expected_node_feature_for_met_sidechain_C3]
        for i, node in enumerate(self.node_features):
            np.testing.assert_array_almost_equal(node, ordered_list_of_expected_node_features[i])
        
        np.testing.assert_array_almost_equal(self.edge_features, expected_bond_features)
        np.testing.assert_array_almost_equal(self.edge_indices, expected_edge_indices)

    def test_partial_charges_heuristic(self):

        #the partial changes should be either positive (1) or negative (-1)
        partial_charges_dict = {"N": -1, # clearly negative since it should pull e from alpha C
                                "alpha_C": 1, # clearly positive since it should donate e to N
                                "met_ketone_C": 1, # clearly positive since it should donate e to O
                                "met_ketone_O": -1} # clearly negative since it should pull e from C
        #assuming that the ordering of the nodes is
        # N, alpha_C, met_ketone_C, met_ketone_O, met_sidechain_C1, met_sidechain_C2, met_sidechain_S, met_sidechain_C3
        for i, node in enumerate(self.node_features):
            if i < 3:
                expected_patial_charge = list(partial_charges_dict.values())[i]
                if expected_patial_charge > 0:
                    self.assertGreater(node[10], 0)
                else:
                    self.assertLess(node[10], 0)
class TestGetProteinOrLigandIds(unittest.TestCase):
    def test_get_protein_or_ligand_ids_methionine(self):
        protein_or_ligand_ids = get_protein_or_ligand_ids(f"{HOME_DIR}/test_data/pdb/test_methionine.pdb")
        #we know that this receptor is just the methionine molecule
        #the first 8 atoms are the protein atoms
        expected_protein_or_ligand_ids = [-1, -1, -1, -1, -1, -1, -1, -1]
        
        np.testing.assert_array_equal(protein_or_ligand_ids, expected_protein_or_ligand_ids)

    def test_get_protein_or_ligand_ids_actual_complex(self):
        protein_or_ligand_ids = get_protein_or_ligand_ids(f"{HOME_DIR}/test_data/pdb/test_101mA_complex.pdb")
        #the first 1221 atoms are protein atoms 
        expected_protein_or_ligand_ids = [-1] * 1221
        #the next 43 atoms are ligand atoms
        expected_protein_or_ligand_ids.extend([1] * 43)

        np.testing.assert_array_equal(protein_or_ligand_ids, expected_protein_or_ligand_ids)

class TestPocketGenerator(unittest.TestCase):
    def test_read_molecule(self):
        protein = read_molecule(f"{HOME_DIR}/test_data/pdb/test_101mA_complex.pdb")
        self.assertEqual(len(protein.atoms), 1264)
    
    def test_find_pocket_atoms(self):
        protein = read_molecule(f"{HOME_DIR}/test_data/pdb/1a0q_protein.pdb")
        ligand = read_molecule(f"{HOME_DIR}/test_data/pdb/1a0q_ligand.mol2")
        pocket_atoms = find_pocket_atoms_OB(protein.atoms, ligand.atoms, 5)
        print(len(pocket_atoms))

    def test_find_and_save_pocket_atoms(self):
        find_and_save_pocket_atoms(f"{HOME_DIR}/test_data/pdb/1a0q_protein.pdb", 
                                   f"{HOME_DIR}/test_data/pdb/1a0q_ligand.mol2", 
                                   5, 
                                   f"{HOME_DIR}/test_data/pdb/1a0q_gen_pocket.pdb")

def main():
    unittest.main()

if __name__ == '__main__':
    main()
