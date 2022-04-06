from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sequence_processing import remove_duplicates

#read in a genome, carry out six frame translation, return all ORFs + write them to a file
def read_genome(genomefile, table, minimum_protein_length, outputfile):

    record = SeqIO.read(genomefile, "fasta")
    
    all_records = []

    for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(record)-frame) // 3) #Multiple of three
            for pro in nuc[frame:frame+length].translate(table).split("*"):
                if len(pro) >= minimum_protein_length:
                    seq_record = SeqRecord(pro)
                    seq_record.id = str(len(pro))
                    seq_record.description = str(record.id) + " " + str(record.description)
                    all_records.append(seq_record)
                    all_records.append(seq_record)
         
    all_records = remove_duplicates(all_records)
    
    SeqIO.write(all_records, outputfile, "fasta")

    return(all_records)

#read in a list of proteins
def read_proteins(proteinfile):

    all_proteins = []

    for protein in SeqIO.parse(proteinfile, "fasta"):
        all_proteins.append(protein) 

    return all_proteins