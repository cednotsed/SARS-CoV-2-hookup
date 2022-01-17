import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO, Seq
from pathlib import Path
from dataclasses import dataclass
from typing import List
cwd = Path.cwd()  # Edit paths as appropriate

# Load reference
ref = SeqIO.read(cwd / 'Wuhan-Hu-1_NC_045512.2.fasta', format='fasta')

# Load genome annotations
meta = pd.read_csv(cwd / 'wuhan-hu-1_genome_annotations_V2.csv')
meta = meta.drop('Unnamed: 8', axis=1)


def get_protein(nuc_start, nuc_end):
    seq = ref.seq[nuc_start - 1:nuc_end]   # nuc_start, nuc_end are one-based indices
    protein_sequence = seq.translate()
    return protein_sequence


# Convert gene annotation positions to per position dataframe
final_df = pd.DataFrame()

# For parsing ORF1ab codons
orf1ab = meta.loc[np.logical_and(meta.region == 'ORF1ab', meta.region_type == 'coding'), 'protein_name']
orf1b = ['NSP12a (part 2)', 'NSP13', 'NSP14', 'NSP15a', 'NSP16']
orf1a = list(set(orf1ab) - set(orf1b))
orf1ab_length = 0
orf1a_length = 0
orf1b_length = 4401

for i in range(len(meta)):
    row = meta.iloc[i:i + 1, :]
    genome_start = row.genome_start.iloc[0]
    genome_end = row.genome_end.iloc[0]
    region_type = row.region_type.iloc[0]
    region = row.region.iloc[0]
    nucleotide_length = row.nucleotide_length.iloc[0]
    protein_name = row.protein_name.iloc[0]
    exploded_element = pd.concat([row] * row.nucleotide_length.iloc[0], axis=0, ignore_index=True)

    # Add nucleotide positions
    exploded_element.insert(loc=3,
                            column='nucleotide_pos',
                            value=np.arange(genome_start, genome_end + 1))

    # Add reference nucleotide
    exploded_element.insert(loc=3,
                            column='ref_nuc',
                            value=ref[genome_start - 1:genome_end])  # convert to zero-based to index sequence

    if region_type == 'coding':
        protein_length = int(row.protein_length.iloc[0])

        # Add position of nucleotide in codon
        exploded_element.insert(loc=3,
                                column='pos_in_codon',
                                value=[1, 2, 3] * protein_length)

        # Add amino acid
        protein_sequence = get_protein(genome_start, genome_end)

        exploded_element.insert(loc=3,
                                column='ref_AA',
                                value=[AA for AA in list(protein_sequence) for _ in range(3)])

    elif region_type == 'non-coding':
        for col in ['pos_in_codon', 'codon_number', 'codon_from_gene_start', 'ref_AA']:
            exploded_element.insert(loc=3,
                                    column=col,
                                    value=-1)
    else:
        raise TypeError('Did you add a new region_type?')

    # Non-continuous codon number from start of each protein
    if region_type == 'coding':
        if protein_name == 'NSP12a (part 2)' and region == 'ORF1ab':
            exploded_element.insert(loc=3,
                                    column='codon_number',
                                    value=[pos for pos in np.arange(protein_length) + 1 + 9 for _ in range(3)])
        else:
            exploded_element.insert(loc=3,
                                    column='codon_number',
                                    value=[pos for pos in np.arange(protein_length) + 1 for _ in range(3)])
    else:
        print(region_type, region, protein_name)

    # Add codon position with ORF1ab in continuous codon sequence accounting for -1 frameshift
    codon_start = 0

    if region_type == 'coding' and region == 'ORF1ab':
        # ORF1a
        if protein_name in orf1a:
            exploded_element.insert(loc=3,
                                    column='codon_from_gene_start',
                                    value=[pos for pos in np.arange(protein_length) + 1 + orf1a_length
                                           for _ in range(3)])

            if protein_name != "NSP11":
                orf1a_length += protein_length

        # Account for -1 frameshift
        elif protein_name in orf1b:
            exploded_element.insert(loc=3,
                                    column='codon_from_gene_start',
                                    value=[pos for pos in np.arange(protein_length) + 1 + orf1b_length
                                           for _ in range(3)])

            orf1b_length += protein_length

        else:
            raise RuntimeError('error in ORF1ab continuous codon number!')

    elif region_type == 'coding' and region != 'ORF1ab':
        exploded_element.insert(loc=3,
                                column='codon_from_gene_start',
                                value=[pos for pos in np.arange(protein_length) + 1
                                       for _ in range(3)])

    # Make dataframes for each possible nuc substitution
    df_list = [exploded_element.copy() for _ in range(4)]
    nuc_list = ['A', 'T', 'G', 'C']

    for df, nuc in zip(df_list, nuc_list):
        df.insert(loc=8, column='var_nuc', value=nuc)
        df.insert(loc=5, column='var_AA', value=-1)

        # Iterate through each codon to get variant AA
        if region_type == 'coding':
            assert len(df) % 3 == 0
            num_codons = len(df) // 3
            idx = 0

            for _ in range(num_codons):
                ori_seq = df.loc[idx:idx + 2, 'ref_nuc'].copy()
                for j in range(3):
                    var_seq = ori_seq.copy().reset_index(drop=True)
                    var_seq.loc[j] = df.loc[idx + j, 'var_nuc']
                    var_codon = Seq.Seq(var_seq.str.cat())
                    var_AA = var_codon.translate()
                    df.loc[idx + j, 'var_AA'] = var_AA

                idx += 3

        exploded_2_element = pd.concat(df_list, axis=0)

    # Append
    final_df = final_df.append(exploded_2_element)


# Filter off entries where ref and var nucleotides are the same
# final_df = final_df.loc[final_df.ref_nuc != final_df.var_nuc, :]

# Add mutation type
final_df.loc[:, 'mutation_type'] = -1
final_df.loc[final_df.loc[:, 'ref_AA'] == final_df.loc[:, 'var_AA'], 'mutation_type'] = 'S'
final_df.loc[final_df.loc[:, 'ref_AA'] != final_df.loc[:, 'var_AA'], 'mutation_type'] = 'NS'

# Save
final_df.to_csv(cwd / 'SARS-CoV-2_hookup_table_V3.csv', index=False)




