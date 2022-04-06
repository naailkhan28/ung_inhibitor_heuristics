from Bio import SeqIO

input_file = "lenienthits.fasta"
output_folder = "FASTAs/"

for seq_record in SeqIO.parse(input_file, "fasta"):
    SeqIO.write(seq_record, output_folder + seq_record.description + ".fasta", "fasta")