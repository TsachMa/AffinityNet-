import os
from chimerax.core.commands import run

#directories = ['/Users/tsachmackey/dfs/affinity_net/PDBbind/v2020-other-PL']
directories = ['test_data/pdb/']
counter = 0 

CHIMERA_CONVERT_TO_MOL2 = True

for directory in directories:
    print(sorted(list(os.listdir(directory))))
    for sub_dir in sorted(list(os.listdir(directory))):
        #check subdir is a directory 
        if os.path.isdir(os.path.join(directory, sub_dir)): 
            for name in os.listdir(os.path.join(directory, sub_dir)):
                name = os.path.join(directory, sub_dir, name)
                if name.endswith("_protein.pdb"):            
                    run(session, f"open {name}")

                    # Add hydrogens
                    run(session, "addh")

                    #delete solvent
                    run(session, "delete solvent")

                    # if CHIMERA_CONVERT_TO_MOL2:
                    #     # Convert to mol2
                    #     protein = name.replace(".pdb", ".mol2")
                    #     run(session, f"save {protein}")
                    # else:
                    #     protein = name.replace(".pdb", "_addH.pdb")
                    #     run(session, f"save {protein}")

                    # # Close all opened models
                    # run(session, "close session")
                    
                    ligand_mol2 =  name.replace("_protein.pdb", "_ligand.mol2")

                    # run(session, f"open {protein}")
                    run(session, f"open {ligand_mol2}")

                    run(session, 'combine #0 #1 #2')

                    models = session.models.list()
                    for idx, model in enumerate(models):
                        print(f"Model #{idx}: {model.name}, Type: {type(model)}")

                    run(session, "close #0,1,2,3")  # Closes models with IDs 1, 2, and 3

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

                    run(session, f"save {complex} #4")

                    # Close all opened models
                    run(session, "close session")

                    counter += 1
                    print(f"Processed {counter} files")
