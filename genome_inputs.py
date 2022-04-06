from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

#read in a list of genomes, carry out six frame translation, return all ORFs + write them to a file
def read_genomes_bulk(genomefile, table, minimum_protein_length, outputfile):

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