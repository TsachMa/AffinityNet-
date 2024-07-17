from multiprocessing import Pool, cpu_count, current_process, Process
import os
from preprocessing.shared_state import load_directories, save_directories

PDB_DIRECTORY = 'test_data/pdb/'
CLEAN = True
PARALLEL = True

def process_subdir(list_of_subdirs):
    #execute /Applications/ChimeraX-1.8.app/Contents/MacOS/ChimeraX --nogui --script preprocessing/add_H_and_mol2_chimerax.py sub_dir
    os.system(f"/Applications/ChimeraX-1.8.app/Contents/MacOS/ChimeraX --nogui --script preprocessing/add_H_and_mol2_chimerax.py")

if __name__ == '__main__':

    list_of_subdirs = [sub_dir for sub_dir in sorted(list(os.listdir(PDB_DIRECTORY))) if os.path.isdir(os.path.join(PDB_DIRECTORY, sub_dir))]
    
    print(f"Directories to process: {list_of_subdirs}")

    # Initialize the shared state with the list of directories
    directories = [{'name': sub_dir, 'processed': False} for sub_dir in list_of_subdirs]
    save_directories(directories)

    if CLEAN:
        for sub_dir in list_of_subdirs:
            for name in os.listdir(os.path.join(PDB_DIRECTORY, sub_dir)):
                name = os.path.join(PDB_DIRECTORY, sub_dir, name)
                if name.endswith("_complex.mol2"):
                    os.remove(name)

    if PARALLEL:
        #set the number of cores 
        num_cpus = cpu_count()
        with Pool(processes=num_cpus//2) as pool: 
            pool.map(process_subdir, list_of_subdirs) #list of subdirectories is jsut to manage the number of processes
            
    else:
        for sub_dir in list_of_subdirs:
            process_subdir(sub_dir)