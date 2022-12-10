# https://www.ncbi.nlm.nih.gov/dbvar/content/tools/entrez/

# NAME: grna_grader.py

"""
AUTHOR: Edward Josh Hermano, Moritz Pistauer

  ============== VARIABLE, FUNCTION, etc. NAMING CONVENTIONS ==================
<ALL CAPITOL LETTERS>:  Indicates a symbol defined by a
        #define statement or a Macro.

   <Capitalized Word>:  Indicates a user defined global var, fun, or typedef.

   <all small letters>:  A variable or built in functions.


========================== MODIFICATION HISTORY ==============================
MM/DD/22:
    MOD:     XXX
    AUTHOR:  Edward Josh Hermano
    COMMENT: No comments

12/09/22:
    MOD:     Implement functions to assess the quality of gRNA sequences
    AUTHOR:  Moritz Pistauer
    COMMENT: No comments

====================== END OF MODIFICATION HISTORY ============================
"""

from Bio import Entrez
from Bio import SeqIO
import xmltodict


def check_for_repetitive_bases(grna):
    base_counts = {"A": 0, "C": 0, "G": 0, "U": 0}

    for base in grna:
        base_counts[base] += 1

    for base, count in base_counts.items():
        if count >= 5:
            return True

    return False


# Function to calculate the GC content of a gRNA
def calculate_gc_content(grna):
    # Count the number of G and C nucleotides in the gRNA
    gc_count = grna.count('G') + grna.count('C')

    # Calculate the GC content as a percentage
    gc_content = gc_count / len(grna)

    return gc_content

# Function to check for off-target sites in a gRNA sequence
def off_target_sites_found(grna, genome):
    # Loop through every substring of the same length as the gRNA in the reference genome
    for i in range(len(genome) - len(grna) + 1):
        substring = genome[i:i + len(grna)]
        # If the gRNA sequence matches the substring, return True
        if grna == substring:
            return True

    # If no matching substring is found, return False
    return False


# Function to assess the quality of a gRNA using the metrics in the CRISPR-CAS9 slide deck
def quality_of_grna(grna):
    # Check the length of the gRNA
    if len(grna) >= 20:
        print('Invalid: gRNA must be 20 or more nucleotides.')
        return False

    # Check the PAM sequence for C's and T's
    pam_sequence = grna[-3:]
    if 'C' in pam_sequence or 'T' in pam_sequence:
        print('Invalid: pam sequnece contains C or T')
        return False

    # Check for off-target sites (not implemented in this example)
    if off_target_sites_found(grna):
        print('Invalid: gRNA has off-target sites.')
        return False

    # Check for GC content
    gc_content = calculate_gc_content(grna)
    if gc_content < 0.3 or gc_content > 0.7:
        print('Invalid: gRNA has GC content outside the acceptable range 30% - 70%.')
        return False

    # Check for repetitive bases like AAAAA, CCCCC, GGGGG, UUUUU
    if check_for_repetitive_bases(grna):
        print('gRNA continues repetitive bases of length 5')
        return False

    return True


def read_grna_file():
    # Read grna file, store only sequences
    seq_reader = open('grna_sequences.txt')
    seq_lines = []

    # NOTE: idk how y'all want the logic on how to format/read the txt file, but this works for the last programming assignment
    
    for line in seq_reader:
        if line[0:4] == '>Seq':
            continue
        elif line[0:1] in ['A', 'T', 'G', 'C', 'N', 'a', 't', 'g', 'c', 'n']:
            seq_lines.append(line.strip().upper())

    seq_reader.close()

    return seq_lines
    

#MAIN Function
def main():
    print('+=== START ===+')

    Entrez.email = 'edwardjosh.hermano@sjsu.edu' # provide your email address

    # (From Mod5_Week14_p2 slides, page 17) (Also on Mod5_Week15_p2)
    # gRNA Design Algorithm
    # Given a DNA Sequence, design efficient gRNAs
    # 1) Find as many gRNAs as possible
    # 2) Map each gRNA to a Reference Genome
    # 3) Score each gRNA


    ### 1) Find as many gRNAs as possible ###
    # Read file to get list of sequences
    print("Fetching sequences of gRNAs from file 'grna_sequences.txt'")
    grna_seqs = read_grna_file()
    

    # Print sequences
    print('##########')
    print('Sequence data:')
    for seq in grna_seqs:
        print(seq)
    print()

    ### 2) Map each gRNA to a Reference Genome ###

    print('Fetching Homo Sapien Chromosome 7...')
    # Download homo sapien genome to cross reference with gRNAs
    # CFTR is located on Chromosome 7 (NC_000007.14) at 117480025..117668665
    
    # https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.ESearch__ESummaryEFetch
    # https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
    
    # nuccore = Entrez Database for Nucleotides
    # NC_000007.14 = NCBI RefSeq for Homo Sapien Chromosome 7

    file_type = 'fasta' # either 'fasta' or 'gb' (GenBank)    ... not sure which we want yet

    handle = Entrez.efetch(db = 'nuccore', id='NC_000007.14', rettype=file_type)
    chromosome = SeqIO.read(handle, format=file_type)
    handle.close()
    # chromosome is of the Seq class from Biopython
    # https://biopython.org/wiki/Seq

    ## Don't print chomosome.seq. It's insanely way too long.
    # print('chromosome.id', chromosome.id)
    # print('chromosome.seq:', chromosome.seq)
    print('chromosome seven:', chromosome)


    # 3) Score each gRNA
    

    print('+===  END  ===+')


# MAKE SURE THIS MODULE IS EXECUTED AND NOT IMPORTED
if __name__ == '__main__':
    main()
else:
    print('Module grna_grader is intended to be executed and not imported.')