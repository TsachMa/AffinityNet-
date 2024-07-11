import os 

#get the directory of the current file
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# get the root directory of the project (2 levels up from the current file)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, os.pardir, os.pardir))

ES_THRESHOLD = 25 #kJ/mol
POCKET_THRESHOLD = 10 #Å

# Van der Waals radii for common elements (in pm)
VDW_RADII = {
    'H': 120, 'He': 140, 'Li': 182, 'Be': 153, 'B': 192, 'C': 170, 'N': 155,
    'O': 152, 'F': 147, 'Ne': 154, 'Na': 227, 'Mg': 173, 'Al': 184, 'Si': 210,
    'P': 180, 'S': 180, 'Cl': 175, 'Ar': 188, 'K': 275, 'Ca': 231, 'Sc': 211,
    'Ti': 200, 'V': 200, 'Cr': 200, 'Mn': 200, 'Fe': 200, 'Co': 200, 'Ni': 163,
    'Cu': 140, 'Zn': 139, 'Ga': 187, 'Ge': 211, 'As': 185, 'Se': 190, 'Br': 185,
    'Kr': 202, 'Rb': 303, 'Sr': 249, 'Y': 200, 'Zr': 200, 'Nb': 200, 'Mo': 200,
    'Tc': 200, 'Ru': 200, 'Rh': 200, 'Pd': 163, 'Ag': 172, 'Cd': 158, 'In': 193,
    'Sn': 217, 'Sb': 206, 'Te': 206, 'I': 198, 'Xe': 216, 'Cs': 343, 'Ba': 268,
}

#NOTE: in the case where there are 2 possible oxidation states, we take the larger one since 
#the predictive penality of incorrectly adding less edges is smaller than the one associated with incorrectly 
# adding more edges

METAL_OX_STATES = {
    'Li': 1,
    'Na': 1,
    'K': 1,
    'Rb': 1,
    'Cs': 1,
    'Be': 2,
    'Mg': 2,
    'Ca': 2,
    'Sr': 2,
    'Ba': 2,
    'Fe': 3, #could be [2, 3],
    'Cu': 2, #could be [1, 2],
    'Zn': 2,
    'Ni': 2,
    'Co': 3, #could be [2, 3],
    'Cr': 6, #could be [3, 6],
    'Mn': 7, #could be [2, 7],
    'Ag': 2, #could be [1, 2],
    'Au': 3,
    'Pt': 4, #could be [2, 4],
    'Hg': 2, #could be [1, 2]
}