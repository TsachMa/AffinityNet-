#!/usr/bin/env bash
export CHIMERA="/Applications/Chimera.app/Contents/Resources"
export PATH="$CHIMERA/bin:$PATH"
export PYTHONPATH="$CHIMERA/lib:$PYTHONPATH"
chimerax --nogui --script preprocessing/add_H_and_mol2_chimerax.py