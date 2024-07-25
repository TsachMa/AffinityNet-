from time import time
import os
from tqdm import tqdm
from rdkit import Chem
from contextlib import contextmanager

@contextmanager
def time_block(label):
    start = time()
    try:
        yield
    finally:
        end = time()
        print(f"{label}: {end - start:.2f} seconds")