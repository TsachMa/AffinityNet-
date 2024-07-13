import os
import sys
from multiprocessing import Pool, cpu_count, Manager
from chimerax.core.commands import run

def process_file(file_path, failed_files):
    global counter
    CHIMERA_CONVERT_TO_MOL2 = True

    try:
        run(session, f"open {file_path}")

        # Add hydrogens
        run(session, "addh")

        # Delete solvent
        run(session, "delete solvent")

        ligand_mol2 = file_path.replace("_protein.pdb", "_ligand.mol2")

        run(session, f"open {ligand_mol2}")

        run(session, 'combine #0 #1')

        models = session.models.list()
        for idx, model in enumerate(models):
            print(f"Model #{idx}: {model.name}, Type: {type(model)}")

        if CHIMERA_CONVERT_TO_MOL2:
            # Convert to mol2
            complex = file_path.replace("_protein.pdb", "_complex.mol2")
        else:
            complex = file_path.replace("_protein.pdb", "_complex.pdb")

        # Calculate AMBER charges for the combined complex
        run(session, "addcharge method am1-bcc")

        models = session.models.list()
        for idx, model in enumerate(models):
            print(f"Model #{idx}: {model.name}, Type: {type(model)}")

        run(session, f"save {complex} #3")

        # Close all opened models
        run(session, "close session")

        counter += 1
        print(f"Processed {counter} files")
        
    except Exception as e:
        print(f"An error occurred while processing {file_path}: {e}")
        failed_files.append(file_path)
        run(session, "close all")

def process_directory(directory):
    files = []
    for sub_dir in sorted(os.listdir(directory)):
        sub_dir_path = os.path.join(directory, sub_dir)
        if os.path.isdir(sub_dir_path):
            for name in os.listdir(sub_dir_path):
                if name.endswith("_protein.pdb"):
                    file_path = os.path.join(sub_dir_path, name)
                    files.append(file_path)
    return files

def process_file_wrapper(args):
    process_file(*args)

if __name__ == "__main__":
    from chimerax.core.session import Session
    session = Session()

    directories = sys.argv[1:]
    all_files = []
    for directory in directories:
        all_files.extend(process_directory(directory))

    num_cores = cpu_count()

    # Use Manager to keep track of failed files
    with Manager() as manager:
        failed_files = manager.list()
        pool = Pool(num_cores)

        # Use a wrapper to pass the shared list to the pool workers
        pool.map(process_file_wrapper, [(file, failed_files) for file in all_files])

        pool.close()
        pool.join()

        # Print the list of failed files
        if failed_files:
            print("The following files were not successfully processed:")
            for file in failed_files:
                print(file)
