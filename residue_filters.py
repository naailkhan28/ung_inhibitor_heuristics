from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sequence_processing import remove_duplicates
import re

#filter for glycine + proline count range
def glyprofilter(proteinlist, mincount, maxcount, outputfile):
    """
    Filter a list of proteins based on the count of glycine (G) and proline (P) residues falling within a specified range.

    The function filters proteins based on the following criteria:
    - The count of glycine and proline residues in each protein falls within the range [mincount, maxcount].

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.
        mincount (int): Minimum count of glycine and proline residues.
        maxcount (int): Maximum count of glycine and proline residues.
        outputfile (str): Path to the output file where filtered proteins will be written in FASTA format.

    Returns:
        list: List of SeqRecord objects representing the proteins that have the count of glycine and proline residues falling within the specified range.

    """

    correct_glypro = []
    
    for protein in proteinlist:
        count = ProteinAnalysis(str(protein.seq)).count_amino_acids()["P"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["G"]

        if count >= mincount and count <= maxcount:
            correct_glypro.append(protein)

    SeqIO.write(correct_glypro, outputfile, "fasta")

    return(correct_glypro)

#filter for min acidic:basic ratio
def chargeratiofilter(proteinlist, minratio, outputfile):
    """
    Filter a list of proteins based on the acidic:basic amino acid ratio.

    The function filters proteins based on the following criteria:
    - The ratio of acidic (E, D, C) amino acids to basic (K, R, H) amino acids is greater than minratio.
    - If the basic amino acid count is zero, the protein is included in the filtered list regardless of the ratio.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.
        minratio (float): Minimum acidic:basic ratio required for inclusion.
        outputfile (str): Path to the output file where filtered proteins will be written in FASTA format.

    Returns:
        list: List of SeqRecord objects representing the proteins that meet the acidic:basic ratio criterion.

    """

    correct_ratio = []

    for protein in proteinlist:
        acidcount = ProteinAnalysis(str(protein.seq)).count_amino_acids()["E"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["D"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["C"]
        basecount = ProteinAnalysis(str(protein.seq)).count_amino_acids()["K"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["R"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["H"]

        try:
            if acidcount / basecount > minratio:
                correct_ratio.append(protein)
        except ZeroDivisionError:
            correct_ratio.append(protein)

    SeqIO.write(correct_ratio, outputfile, "fasta")

    return(correct_ratio)

#filter for longest acidic-residue-free region
def acidfree(proteinlist, minstretch, outputfile):
    """
    Filter a list of proteins based on the length of the longest acidic-residue-free region.

    The function filters proteins based on the following criteria:
    - The length of the longest stretch of consecutive amino acids without any acidic residues (E, D) is greater than or equal to minstretch.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.
        minstretch (int): Minimum length of the longest acidic-residue-free region.
        outputfile (str): Path to the output file where filtered proteins will be written in FASTA format.

    Returns:
        list: List of SeqRecord objects representing the proteins that have a sufficiently long acidic-residue-free region.

    """

    correct_acidfree = []

    acidic = ["E", "D"]

    for protein in proteinlist:

        sequence = str(protein.seq)

        for i in range(0, len(sequence) - minstretch):
            window = sequence[i:i+minstretch]

            has_acidic = False

            for x in window:
                if x in acidic: 
                    has_acidic = True

        if has_acidic:
            correct_acidfree.append(protein)
    
    correct_acidfree = remove_duplicates(correct_acidfree)

    SeqIO.write(correct_acidfree, outputfile, "fasta")

    return(correct_acidfree)

#filter for number of ED motifs
def edfinder(proteinlist, minmotifs, outputfile):
    """
    Filter a list of proteins based on the count of ED motifs.

    The function filters proteins based on the following criteria:
    - The count of ED motifs (EE, DD, DE, ED) in each protein is greater than or equal to minmotifs.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.
        minmotifs (int): Minimum count of ED motifs.
        outputfile (str): Path to the output file where filtered proteins will be written in FASTA format.

    Returns:
        list: List of SeqRecord objects representing the proteins that have the count of ED motifs greater than or equal to minmotifs.

    """

    correct_motifs = []

    for protein in proteinlist:

        z_sequence = re.sub("(EE|DD|DE|ED)", "Z", str(protein.seq))

        ed_list = re.findall("Z", z_sequence)
        
        if len(ed_list) >= 3:
            correct_motifs.append(protein)

    SeqIO.write(correct_motifs, outputfile, "fasta")

    return(correct_motifs)