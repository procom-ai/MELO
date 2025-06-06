from muliti_chain_score import *
import os
from Bio.PDB import PDBParser, PDBIO
import pandas as pd
pdb_name1 = f'6A6M'
pdb_name2 = f'6A6N'
# pdb_file1 = os.path.join('/mnt/h/MELO/NC/low_sequence_similar', f"{pdb_name1}.pdb")
# pdb_file2 = os.path.join('/mnt/h/MELO/NC/low_sequence_similar', f"{pdb_name2}.pdb")
# dssp_file1 = os.path.join('/mnt/h/MELO/NC/low_sequence_similar', f"{pdb_name1}.dssp")
# dssp_file2 = os.path.join('/mnt/h/MELO/NC/low_sequence_similar', f"{pdb_name2}.dssp")
# chainid1=['A','C']
# chainid2=['A','B']
# aligned_seq1, aligned_seq2, aligned_seq1_ss, aligned_seq2_ss, aligned_seq1_vector, aligned_seq2_vector, score_engry, score_engry_count, mean_score, mean_score_rmsd, d_rmsd = aligh_match(
#     pdb_file1, pdb_name1, dssp_file1, pdb_file2, pdb_name2, dssp_file2, chainid1, chainid2 )
# aligned_seq1_list = list(aligned_seq1)
# aligned_seq2_list = list(aligned_seq2)
# d_rmsd = pd.DataFrame(d_rmsd)
# assert len(aligned_seq1_list) == len(d_rmsd)
# assert len(aligned_seq2_list) == len(d_rmsd)
# d_rmsd.insert(0, 'Aligned_Seq1', aligned_seq1_list)
# d_rmsd.insert(1, 'Aligned_Seq2', aligned_seq2_list)
# d_rmsd.to_csv('/mnt/h/MELO/NC/review2/residue/DRMSD/6A6M_6A6N.csv')
d_rmsd = pd.read_csv('/mnt/h/MELO/NC/review2/residue/DRMSD/6A6M_6A6N.csv')
aligned_seq1_list = d_rmsd['Aligned_Seq1'].tolist()
aligned_seq2_list = d_rmsd['Aligned_Seq2'].tolist()
def color_structure_by_mean_rmsd(pdb_file, mean_rmsd_list, aligned_seq, chain_ids, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)

    i = 0  
    for model in structure:
        for chain in model:
            if chain.id in chain_ids:
                for residue in chain:
                    while i < len(mean_rmsd_list) and aligned_seq[i] == '-':
                        i += 1
                    if i >= len(mean_rmsd_list):
                        break
                    value = mean_rmsd_list[i]
                    if value is not None:
                        for atom in residue:
                            atom.set_bfactor(value * 100)
                    i += 1

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)


d_rmsd_values = d_rmsd.iloc[:, 2:]


j = 200
k = 280


mean_rmsd = d_rmsd_values.iloc[:, j:k+1].mean(axis=1)


d_rmsd["mean_rmsd"] = mean_rmsd

normalized_rmsd = pd.Series(mean_rmsd) / pd.Series(mean_rmsd).max()
normalized_rmsd = normalized_rmsd.fillna(0).tolist()

color_structure_by_mean_rmsd(
    pdb_file="/mnt/h/MELO/NC/low_sequence_similar/6A6N.pdb",
    mean_rmsd_list=normalized_rmsd,
    aligned_seq=aligned_seq2_list,  
    chain_ids=["A", "B"],
    output_pdb="/mnt/h/MELO/NC/review2/residue/DRMSD/PDB/6A6N_colored_A2.pdb"
)
color_structure_by_mean_rmsd(
    pdb_file="/mnt/h/MELO/NC/low_sequence_similar/6A6M_ZZ.pdb",
    mean_rmsd_list=normalized_rmsd,
    aligned_seq=aligned_seq1_list,  
    chain_ids=["A", "C"],
    output_pdb="/mnt/h/MELO/NC/review2/residue/DRMSD/PDB/6A6M_colored_A2.pdb"
)

#pymol
#set_color p2,[70, 70, 196]   set_color p2,[51, 51, 204] 
# set_color no,[200,200,200]
# set_color p1,[200, 0, 0]
# set_color target,[243, 207, 71]
# spectrum b, no_p2, 6A6M_colored_A1, minimum=0, maximum=100
# spectrum b, no_p1, 6A6N_colored_A1, minimum=0, maximum=100
# color target,sele 
# set cartoon_transparency,0.8,sele
