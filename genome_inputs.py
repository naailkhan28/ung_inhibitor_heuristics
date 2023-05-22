from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

#read in a list of genomes, carry out six frame translation, return all ORFs + write them to a file
def read_genomes_bulk(genomefile, table, minimum_protein_length, outputfile):
    """
    Read in a list of genomes, perform six-frame translation, and return all Open Reading Frames (ORFs).

    Args:
        genomefile (str): Path to the file containing the genomes in FASTA format.
        table (int): Integer representing the translation table to use.
        minimum_protein_length (int): Minimum length of proteins to consider as ORFs.
        outputfile (str): Path to the output file where the ORFs will be written in FASTA format.

    Returns:
        list: List of SeqRecord objects representing the ORFs found.

    """

    all_records = []

    for record in SeqIO.parse(genomefile, "fasta"):
        for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
            for frame in range(3):
                length = 3 * ((len(record)-frame) // 3) #Multiple of three
                for pro in nuc[frame:frame+length].translate(table).split("*"):
                    if len(pro) >= minimum_protein_length:
                        seq_record = SeqRecord(pro)
                        seq_record.id = str(len(pro))
                        seq_record.description = str(record.id) + " " + str(record.description)
                        all_records.append(seq_record)
         
    SeqIO.write(all_records, outputfile, "fasta")

    return(all_records)

#split a list of genomes into two lists based on size
def split_genomes_size(genomefile, cutoff, highoutput, lowoutput):
    """
    Split a list of genomes into two lists based on their size.

    Genomes larger than the specified cutoff will be added to the 'big_genomes' list,
    while genomes smaller than or equal to the cutoff will be added to the 'small_genomes' list.

    Args:
        genomefile (str): Path to the file containing the genomes in FASTA format.
        cutoff (int): Size threshold for splitting the genomes into big and small categories.
        highoutput (str): Path to the output file where the big genomes will be written in FASTA format.
        lowoutput (str): Path to the output file where the small genomes will be written in FASTA format.

    Returns:
        None
    """

    big_genomes = []
    small_genomes = []

    for record in SeqIO.parse(genomefile, "fasta"):
        if len(str(record.seq)) > cutoff:
            big_genomes.append(record)
        else:
            small_genomes.append(record)
    
    SeqIO.write(big_genomes, highoutput, "fasta")
    SeqIO.write(small_genomes, lowoutput, "fasta")

#split a list of genomes into three lists based on GC content
def split_genomes_gc(genomefile, highcutoff, lowcutoff, medcutoff1, medcutoff2, medcutoff3, highoutput, medoutput1, medoutput2, medoutput3, medoutput4, lowoutput):
    """
    Split a list of genomes into multiple lists based on their GC content.

    Genomes will be categorized into six lists based on their GC content ranges:
    - 'high_gc' for genomes with GC content greater than 'highcutoff'
    - 'low_gc' for genomes with GC content less than 'lowcutoff'
    - 'med_gc1' to 'med_gc4' for genomes with GC content falling within four intermediate ranges:
        - 'med_gc1' for GC content between 'lowcutoff' and 'medcutoff1'
        - 'med_gc2' for GC content between 'medcutoff1' and 'medcutoff2'
        - 'med_gc3' for GC content between 'medcutoff2' and 'medcutoff3'
        - 'med_gc4' for GC content between 'medcutoff3' and 'highcutoff'

    Args:
        genomefile (str): Path to the file containing the genomes in FASTA format.
        highcutoff (float): Upper threshold for the high GC content range.
        lowcutoff (float): Lower threshold for the low GC content range.
        medcutoff1 (float): Lower threshold for the first intermediate GC content range.
        medcutoff2 (float): Upper threshold for the first intermediate GC content range and lower threshold for the second range.
        medcutoff3 (float): Upper threshold for the second intermediate GC content range and lower threshold for the third range.
        highoutput (str): Path to the output file where genomes with high GC content will be written in FASTA format.
        medoutput1 (str): Path to the output file where genomes in the first intermediate GC content range will be written in FASTA format.
        medoutput2 (str): Path to the output file where genomes in the second intermediate GC content range will be written in FASTA format.
        medoutput3 (str): Path to the output file where genomes in the third intermediate GC content range will be written in FASTA format.
        medoutput4 (str): Path to the output file where genomes in the fourth intermediate GC content range will be written in FASTA format.
        lowoutput (str): Path to the output file where genomes with low GC content will be written in FASTA format.

    Returns:
        None
    """

    high_gc = []
    low_gc = []
    med_gc1 = []
    med_gc2 = []
    med_gc3 = []
    med_gc4 = []

    for record in SeqIO.parse(genomefile, "fasta"):

        percentage = GC(str(record.seq))

        if percentage < lowcutoff:
            low_gc.append(record)
        elif percentage > highcutoff:
            high_gc.append(record)
        elif percentage >= lowcutoff and percentage <= medcutoff1:
            med_gc1.append(record)
        elif percentage >= medcutoff1 and percentage <= medcutoff2:
            med_gc2.append(record)
        elif percentage >= medcutoff2 and percentage <= medcutoff3:
            med_gc3.append(record)
        elif percentage >= medcutoff3 and percentage <= highcutoff:
            med_gc4.append(record)

    SeqIO.write(high_gc, highoutput, "fasta")
    SeqIO.write(med_gc1, medoutput1, "fasta")
    SeqIO.write(med_gc2, medoutput2, "fasta")
    SeqIO.write(med_gc3, medoutput3, "fasta")
    SeqIO.write(med_gc4, medoutput4, "fasta")
    SeqIO.write(low_gc, lowoutput, "fasta")