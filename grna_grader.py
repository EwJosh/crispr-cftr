# https://www.ncbi.nlm.nih.gov/dbvar/content/tools/entrez/

from Bio import Entrez
from Bio import SeqIO
import xmltodict

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