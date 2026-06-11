''' Run all the files inside the specified folders'''
import os
import sys

if __name__ == "__main__":
    # Define the folders containing the scripts to run
    path = os.path.dirname(__file__)

    folders = ["L2_mis", "geometric_mis"]
    
    # Skip files name:
    names_to_skip = [  ]

    folders_paths = [os.path.join(path, folder) for folder in folders]
    for folder in folders_paths:
        for filename in os.listdir(folder):
            if filename.endswith(".py") and filename not in names_to_skip:
                file_path = os.path.join(folder, filename)
                print(f"Running {file_path}...")
                os.system(f"python {file_path}")
    print("All scripts have been executed.")