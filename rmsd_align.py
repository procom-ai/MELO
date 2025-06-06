from Bio.PDB import PDBParser, PDBIO
import numpy as np
import pandas as pd
from Bio.PDB.Polypeptide import is_aa
# ========== pipeline ==========
#/mnt/e/topic/mutant_protein_structure/TMalign 1KH0.pdb 1PGB.pdb > 1KH0_1PGB.txt
#PYMOL  STEP 1: cmd.show_as("cartoon"   ,"all") 
# STEP 2:set_color p1,[200, 0, 0]
# spectrum b, white_p1, 6A6N_colored_A1, minimum=0, maximum=100
#set_color p2,[51, 51, 204] 
# spectrum b, white_p2, 6A6M_colored_A1, minimum=0, maximum=100
pdb_file1 = '/mnt/h/MELO/NC/review2/aligned_pdb_files/aligned_pdb_files/align/6A6M_align.pdb'
pdb_file2 = '/mnt/h/MELO/NC/review2/aligned_pdb_files/aligned_pdb_files/align/6A6N_align.pdb'
tm_output = "/mnt/h/MELO/NC/low_sequence_similar/6A6N_6A6M.txt"

def extract_aligned_sequences_blocks(tm_align_output_file):
    with open(tm_align_output_file, 'r') as f:
        lines = f.readlines()

    blocks = []
    i = 0
    while i < len(lines) - 2:
        line1 = lines[i].strip()
        line2 = lines[i+1].strip()
        line3 = lines[i+2].strip()
        if (
            all(c.isalpha() or c == '-' for c in line1) and
            all(c in ":. " for c in line2) and
            all(c.isalpha() or c == '-' for c in line3)
        ):
            blocks.append((line1, line3, line2))
            i += 3
        else:
            i += 1
    return blocks  

def parse_structure_ca_coords(pdb_file, chain_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    chain = structure[0][chain_id]
    coords = []
    for res in chain:
        if 'CA' in res:
            coords.append((chain.id, res.get_id()[1], res.get_resname(), res['CA'].get_coord()))
    return coords

def compute_rmsd(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

def map_aligned_residues(seq1, seq2, coords1, coords2):
    assert len(seq1) == len(seq2)
    i1 = i2 = 0
    result = []
    for a1, a2 in zip(seq1, seq2):
        if a1 != '-' and a2 != '-':
            chain1, resi1, resn1, coord1 = coords1[i1]
            chain2, resi2, resn2, coord2 = coords2[i2]
            rmsd = compute_rmsd(coord1, coord2)
            result.append({
                "Chain1": chain1, "Resi1": resi1, "Resn1": resn1,
                "Chain2": chain2, "Resi2": resi2, "Resn2": resn2,
                "RMSD": round(rmsd, 3)
            })
        if a1 != '-':
            i1 += 1
        if a2 != '-':
            i2 += 1
    return result

def color_structure_by_rmsd(pdb_file, rmsd_df, chain_ids, output_pdb, proteinid):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    rmsd_dict = {
        (row[f"Chain{proteinid}"], row[f"Resi{proteinid}"]): row["RMSD"]
        for _, row in rmsd_df.iterrows()
    }
    max_rmsd = rmsd_df["RMSD"].max()
    for model in structure:
        for chain in model:
            if chain.id in chain_ids:
                for residue in chain:
                    if not is_aa(residue, standard=True):
                        continue
                    if 'CA' in residue:
                        resi = residue.get_id()[1]
                        key = (chain.id, resi)
                        # print(key)
                        if key in rmsd_dict:
                            norm_rmsd = rmsd_dict[key] / max_rmsd
                            for atom in residue:
                                atom.set_bfactor(norm_rmsd * 100.0)
                        else:
                            for atom in residue:
                                atom.set_bfactor(100.0)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

blocks = extract_aligned_sequences_blocks(tm_output)


seq1_A, seq2_A, _ = blocks[0]
# seq1_B, seq2_B, _ = blocks[0]

coords1_A = parse_structure_ca_coords(pdb_file1, "A")
# coords1_B = parse_structure_ca_coords(pdb_file1, "B")
coords2_A = parse_structure_ca_coords(pdb_file2, "A")
# coords2_B = parse_structure_ca_coords(pdb_file2, "B")

rmsd_A = map_aligned_residues(seq1_A, seq2_A, coords1_A, coords2_A)
# rmsd_B = map_aligned_residues(seq1_B, seq2_B, coords1_B, coords2_B)
# print(seq1_A)

df = pd.DataFrame(rmsd_A)
# df = pd.DataFrame(rmsd_A + rmsd_B)
# print(df)
df.to_csv("/mnt/h/MELO/NC/review2/residue/7OXW_6HKR.csv", index=False)


color_structure_by_rmsd(
    pdb_file=pdb_file2,
    rmsd_df=df,
    chain_ids=["A"],
    output_pdb="/mnt/h/MELO/NC/review2/residue/6HKR_color.pdb",
    proteinid = 2
)

color_structure_by_rmsd(
    pdb_file=pdb_file1,
    rmsd_df=df,
    chain_ids=["A"],
    output_pdb="/mnt/h/MELO/NC/review2/residue/6HZM_color.pdb",
    proteinid = 1
)
