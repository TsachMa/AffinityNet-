import os
import json
from multiprocessing import Lock

PDB_DIRECTORY = 'test_data/pdb/'
DIRECTORY_FILE = os.path.join(PDB_DIRECTORY,'directories.json')
LOCK = Lock()

def load_directories():
    if not os.path.exists(DIRECTORY_FILE):
        return []

    with LOCK:
        with open(DIRECTORY_FILE, 'r') as f:
            return json.load(f)

def save_directories(directories):
    with LOCK:
        with open(DIRECTORY_FILE, 'w') as f:
            json.dump(directories, f)

def get_next_directory():
    with LOCK:
        directories = load_directories()
        for directory in directories:
            if not directory.get('processed'):
                directory['processed'] = True
                save_directories(directories)
                return directory['name']
    return None
