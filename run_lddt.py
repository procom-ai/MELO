from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Polypeptide import is_aa
import pandas as pd
LDDT = pd.read_csv("/mnt/h/MELO/NC/review2/residue/DRMSD/ABCG2_lddt.csv")
lddt_v = LDDT['lddt'].tolist()

def color_structure_by_mean_rmsd(pdb_file, mean_rmsd_list, aligned_seq, chain_ids, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)

    i = 0  
    for model in structure:
        for chain in model:
            if chain.id in chain_ids:
                for residue in chain:
                    if not is_aa(residue):
                         continue
                    # print(residue)
                    while i < len(mean_rmsd_list) and pd.isna(aligned_seq[i]):
                        i += 1
                        # print(i)
                    if i >= len(mean_rmsd_list):
                        break
                    # 设置 bfactor
                    value = mean_rmsd_list[i]
                    print(i,residue,value,len(mean_rmsd_list))
                    if value is not None:
                        for atom in residue:
                            atom.set_bfactor(value * 100)
                            # print(atom)
                    i += 1

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
aligned_seq1_list = LDDT['mdl_res'].tolist()
aligned_seq2_list = LDDT['ref_res'].tolist()
color_structure_by_mean_rmsd(
    pdb_file="/mnt/h/MELO/NC/low_sequence_similar/6HZM.pdb",
    mean_rmsd_list=lddt_v,
    aligned_seq=aligned_seq1_list,  
    chain_ids=["A"],
    output_pdb="/mnt/h/MELO/NC/review2/residue/DRMSD/PDB/6HZM_colored_LDDT.pdb",
)
color_structure_by_mean_rmsd(
    pdb_file="/mnt/h/MELO/NC/low_sequence_similar/6VXH.pdb",
    mean_rmsd_list=lddt_v,
    aligned_seq=aligned_seq2_list, 
    chain_ids=["A"],
    output_pdb="/mnt/h/MELO/NC/review2/residue/DRMSD/PDB/6VXH_colored_LDDT.pdb",
)
