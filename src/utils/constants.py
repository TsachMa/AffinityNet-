import os 

#get the directory of the current file
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# get the root directory of the project (2 levels up from the current file)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, os.pardir, os.pardir))

ES_THRESHOLD = 50 #kJ/mol
POCKET_THRESHOLD = 10 #Å