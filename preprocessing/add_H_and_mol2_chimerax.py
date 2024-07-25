import os
from chimerax.core.commands import run
import argparse
import sys

import os

PDB_DIRECTORY = 'test_data/pdb/'
CHIMERA_CONVERT_TO_MOL2 = True
CHECKPOINTING = True 
LOG_ERRORS = True
PARALLEL = True 

if LOG_ERRORS:
    error_log = open("error_log.txt", "w")

counter = 0

for sub_dir in os.listdir(PDB_DIRECTORY):
    if os.path.isdir(os.path.join(PDB_DIRECTORY, sub_dir)):
        for name in os.listdir(os.path.join(PDB_DIRECTORY, sub_dir)):
            
            name = os.path.join(PDB_DIRECTORY, sub_dir, name)
            
            if name.endswith("_protein.pdb"):

                try: 
                    if CHECKPOINTING:
                        if os.path.exists(name.replace("_protein.pdb", "_complex.mol2")):
                            print(f"Skipping {name}")
                            continue

                    run(session, f"open {name}")

                    # Add hydrogens
                    run(session, "addh")

                    #delete solvent
                    run(session, "delete solvent")
                    
                    ligand_mol2 =  name.replace("_protein.pdb", "_ligand.mol2")

                    # run(session, f"open {protein}")
                    run(session, f"open {ligand_mol2}")

                    run(session, 'combine #0 #1 #2')

                    models = session.models.list()
                    for idx, model in enumerate(models):
                        print(f"Model #{idx}: {model.name}, Type: {type(model)}")

                    if CHIMERA_CONVERT_TO_MOL2:
                        # Convert to mol2
                        complex = name.replace("_protein.pdb", "_complex.mol2")
                    else:
                        complex = name.replace("_protein.pdb", "_complex.pdb")

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
                    print(f"Error processing {name}: {e}")
                    error_log.write(f"Error processing {name}: {e}\n")

run(session, "quit")