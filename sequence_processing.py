from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#turn an input string into a dictionary with index as key and residue as value
def dictify(in_seq):
    """
    Convert an input string into a dictionary with index as the key and residue as the value.

    Args:
        in_seq (str): Input sequence string.

    Returns:
        dict: Dictionary with index as the key and residue as the value.

    """
    
    dict = {}
    count = 0
    for i in in_seq:
        dict[count] = i
        count += 1
    return dict

#replace illegal amino acid residues in protein sequences
def remove_illegal_chars(text):
    """
    Replace illegal amino acid residues in protein sequences.

    Args:
        text (str): Input sequence string.

    Returns:
        str: Sequence string with illegal residues replaced.

    """
    text = text.replace("B", "N")
    text = text.replace("Z", "Q")
    text = text.replace("J", "I")
    text = text.replace("X", "I")
    return text

#remove duplicate sequences
def remove_duplicates(proteinlist):
    """
    Remove duplicate sequences from a list of protein records.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of SeqRecord objects with duplicate sequences removed.

    """
    
    seen = []
    all_records = []

    for record in proteinlist:
        if str(record.seq) not in seen:
            seen.append(str(record.seq))
            all_records.append(record)

    return all_records

#trim all residues before the first putative start codon and remove dupes - takes an input list of potential start codon residues
def trim_residues(proteinlist, startcodons, outputfile):
    """
    Trim all residues before the first putative start codon and remove duplicates.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.
        startcodons (list): List of potential start codon residues.
        outputfile (str): Output file name for writing the trimmed sequences.

    Returns:
        list: List of trimmed SeqRecord objects representing the proteins.

    """

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