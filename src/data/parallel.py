import time
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

def parallel_process_files(directory: str, files: list, process_file, cpu_fraction: float = 0.5):
    """
    Args: 

        directory (str): The path to the directory where the files are located.
        files (list): The list of files to process. Also supports inputs of tuples of files if the function can take that as input 
        process_file (function): The function to apply to each file.

    Description:
        This function finds all files in the given directory that match the given regex pattern and applies the given function to each file concurrently using multiprocessing. 
    """
    
    total_files = len(files)
    converted_files = 0
    start_time = time.time()

    # Use multiprocessing Pool to process files concurrently
    cpus = int(cpu_count() * cpu_fraction)
    with Pool(cpus) as pool:
        for _ in tqdm(pool.imap_unordered(process_file, files), total=total_files):
            converted_files += 1

    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time} seconds")