from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#turn an input string into a dictionary with index as key and residue as value
def dictify(in_seq):
    
    dict = {}
    count = 0
    for i in in_seq:
        dict[count] = i
        count += 1
    return dict

#replace illegal amino acid residues in protein sequences
def remove_illegal_chars(text):
    text = text.replace("B", "N")
    text = text.replace("Z", "Q")
    text = text.replace("J", "I")
    text = text.replace("X", "I")
    return text

#remove duplicate sequences
def remove_duplicates(proteinlist):
    
    seen = []
    all_records = []

    for record in proteinlist:
        if str(record.seq) not in seen:
            seen.append(str(record.seq))
            all_records.append(record)

    return all_records

#trim all residues before the first putative start codon and remove dupes - takes an input list of potential start codon residues
def trim_residues(proteinlist, startcodons, outputfile):

    trimmed_seqs = []

    for protein in proteinlist:

        sequence = remove_illegal_chars(str(protein.seq))
        for i in sequence:
            if i in startcodons:
                pos = sequence.index(i)
                trimmed = sequence[pos:]
                seq_record = SeqRecord(Seq(trimmed))
                seq_record.id = str(protein.id)
                seq_record.description = str(protein.id) + " " + str(protein.description)
                trimmed_seqs.append(seq_record)
                break
    
    SeqIO.write(trimmed_seqs, outputfile, "fasta")

    return(trimmed_seqs)