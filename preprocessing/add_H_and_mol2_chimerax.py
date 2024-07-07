import os
from chimerax.core.commands import run

# Change to folder with data files
#os.chdir("/Users/tsachmackey/dfs/affinity_net/PDBbind/v2020-other-PL")

directories = ['test_data/pdb']

for directory in directories:
    print( sorted(list(os.listdir(directory))))
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

                    # Convert to mol2
                    protein_mol2 = name.replace(".pdb", ".mol2")
                    run(session, f"save {protein_mol2}")

                    # Close all opened models
                    run(session, "close session")
                    
                    ligand_mol2 =  name.replace("_protein.pdb", "_ligand.mol2")

                    run(session, f"open {protein_mol2}")
                    run(session, f"open {ligand_mol2}")

                    run(session, 'combine #0 #1 #2')

                    models = session.models.list()
                    for idx, model in enumerate(models):
                        print(f"Model #{idx}: {model.name}, Type: {type(model)}")

                    # Convert to mol2
                    complex_mol2 = name.replace("_protein.pdb", "_complex.mol2")

                    run(session, f"save {complex_mol2} #3 #4")

                    # Close all opened models
                    run(session, "close session")
