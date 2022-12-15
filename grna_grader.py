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

12/14/22:
    MOD:     Calculate scores
    AUTHOR:  Moritz Pistauer
    COMMENT: No comments

====================== END OF MODIFICATION HISTORY ============================
"""

from Bio import Entrez
from Bio import SeqIO


# Function to assess the quality of a gRNA using the metrics in the CRISPR-CAS9 slide deck
def quality_of_grna(grna):
    # Check the length of the gRNA
    if len(grna) < 20:
        print('gRNA must be 20 or more nucleotides.')
        return False

    # Check the PAM sequence for C's and T's
    pam_sequence = grna[-3:]
    if 'C' in pam_sequence or 'T' in pam_sequence:
        print('pam sequnece contains C or T')
        return False

    # Check for GC content
    gc_count = grna.count('G') + grna.count('C')
    gc_content_percent = gc_count / len(grna)
    if gc_content_percent < 0.3 or gc_content_percent > 0.7:
        print('gRNA has GC content outside the acceptable range 30% - 70%.')
        return False

    return True


def read_grna_file():
    # Read grna file, store only sequences
    seq_reader = open('grna_sequences.txt')
    seq_lines = []

    for line in seq_reader:
        if line[0:4] == '>Seq':
            continue
        elif line[0:1] in ['A', 'T', 'G', 'C', 'N', 'a', 't', 'g', 'c', 'n']:
            seq_lines.append(line.strip().upper())

    seq_reader.close()

    return seq_lines


def count_repetitive_bases(grna):
    base_counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    prev_base = ''
    count = 1

    for base in grna:
        if base == prev_base:
            count += 1
            if count > base_counts[base]:
                base_counts[base] = count
        else:
            count = 1

        prev_base = base

    return base_counts


def score_repetitive_bases(grna_seqs):
    repetition_scores = {}
    for grna in grna_seqs:
        max_rep_count = count_repetitive_bases(grna)
        max_rep = max(max_rep_count.values())
        repetition_scores.setdefault(grna, round(10/max_rep, 2))
    return repetition_scores


# Function to score the GC content of a gRNA
# The lower the GC count is compared to 0.57, the bigger the absolute distance the higher the score
def calculate_gc_content(grna_seqs):
    gc_content_scores = {}
    for grna in grna_seqs:
        gc_count = grna.count('G') + grna.count('C')
        gc_content_percent = gc_count / len(grna)
        gc_content_scores.setdefault(grna, gc_content_percent)

    # Normalize the scores (max_score = 10)
    for gc_score in gc_content_scores:
        distance = abs(gc_content_scores[gc_score]-0.57)/0.57
        gc_content_scores[gc_score] = round(10 * (1 - distance), 2)

    return gc_content_scores


# Position-specific sequence
# Count the number of 'C' and 'T' on the last 10 positions and score it (the more C&T the lower the score)
def score_position(grna_seqs):
    position_scores = {}

    for grna in grna_seqs:
        last_10 = grna[-10:]
        c_count = last_10.count('C')
        t_count = last_10.count('T')
        position_scores.setdefault(grna, c_count + t_count)

    min_count = min(position_scores.values())

    # Normalize the scores (max_score = 10)
    for element in position_scores:
        if position_scores[element] != 0:
            position_scores[element] = round(10*(min_count/position_scores[element]), 2)

    return position_scores


OUTPUT_FILE = 'grna_scores.txt'
# Calculate the total scores
# The results will be written into grna_scores.txt file
def score(grna_seqs):
    overall_scores = {}
    position_scores = score_position(grna_seqs)
    gc_content_scores = calculate_gc_content(grna_seqs)
    repetitive_base_scores = score_repetitive_bases(grna_seqs)

    print('Outputing results file...')
    with open(OUTPUT_FILE, 'w') as f:
        for grna in grna_seqs:
            overall_scores[grna] = round(position_scores[grna] + gc_content_scores[grna] + repetitive_base_scores[grna], 2)
            output = f'\n-------------------------------------' \
                     f'\nsRNA Sequence: {grna}\n\n' \
                     f'Scoring Metrics:\n\tPosition Score: {position_scores[grna]}/10.0\n\t' \
                     f'GC Content Score: {gc_content_scores[grna]}/10.0\n\t' \
                     f'Repetitive Base Score : {repetitive_base_scores[grna]}/10.0\n' \
                     f'==> Assessment Score: {overall_scores[grna]}/30.0\n\n'
            f.write(output)
    print(f'Results output at: {OUTPUT_FILE}')
    return overall_scores
    

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

    for grna in grna_seqs:
        if not quality_of_grna(grna):
            print("gRNA has bad quality: " + grna)

    score(grna_seqs)

    print('+===  END  ===+')


# MAKE SURE THIS MODULE IS EXECUTED AND NOT IMPORTED
if __name__ == '__main__':
    main()
else:
    print('Module grna_grader is intended to be executed and not imported.')