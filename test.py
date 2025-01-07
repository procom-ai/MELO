from melo import melo
protein_file = "TEST.csv"  # Input protein pairs file
pdb_file = "test_pdb"       # Folder containing PDB files
save_folder = "result"         # Folder to save the results
threads = 1                         # Set number of threads
melo(protein_file, pdb_file, save_folder, threads)
