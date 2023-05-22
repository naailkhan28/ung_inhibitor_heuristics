from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#filter below max length and molecular weight
def lengthfilter(proteinlist, maxlength, maxmW, outputfile):
    """
    Filters a list of SeqRecord objects based on the length and molecular weight of their protein sequence.

    Args:
        proteinlist (list): A list of SeqRecord objects.
        maxlength (int): The maximum length of the protein sequence to include.
        maxmW (int): The maximum molecular weight in kDA of the protein sequence to include.
        outputfile (str): The path to the output FASTA file.

    Returns:
        list: A list of filtered SeqRecord objects.

    Raises:
        IOError: If the output file cannot be opened or written.
    """

    correct_size = []

    for protein in proteinlist:
        if ProteinAnalysis(str(protein.seq)).molecular_weight() < maxmW and len(str(protein.seq)) < maxlength:
            correct_size.append(protein)

    SeqIO.write(correct_size, outputfile, "fasta")

    return(correct_size)

#filter for correct pI range
def pifilter(proteinlist, minpi, maxpi, outputfile):
    """
    Filters a list of SeqRecord objects based on the isoelectric point of their protein sequence.

    Args:
        proteinlist (list): A list of SeqRecord objects.
        minpi (float): The minimum isoelectric point of the protein sequence to include.
        maxpi (float): The maximum isoelectric point of the protein sequence to include.
        outputfile (str): The path to the output FASTA file.

    Returns:
        list: A list of filtered SeqRecord objects.

    Raises:
        IOError: If the output file cannot be opened or written.
    """

    correct_pi = []

    for protein in proteinlist:
        if ProteinAnalysis(str(protein.seq)).isoelectric_point() <= maxpi and ProteinAnalysis(str(protein.seq)).isoelectric_point() >= minpi:
            correct_pi.append(protein)

    SeqIO.write(correct_pi, outputfile, "fasta")

    return(correct_pi)

#filter for hydrophobic residue percentage range
def hydrophobicfilter(proteinlist, minvalue, maxvalue, outputfile):
    """
    Filters a list of SeqRecord objects based on the hydrophobicity of their protein sequence.
    Hydrophobicity is defined as the percentage of hydrophobic residues [IVLFMAWT]

    Args:
        proteinlist (list): A list of SeqRecord objects.
        minvalue (float): The minimum hydrophobicity of the protein sequence to include.
        maxvalue (float): The maximum hydrophobicity of the protein sequence to include.
        outputfile (str): The path to the output FASTA file.

    Returns:
        list: A list of filtered SeqRecord objects.

    Raises:
        IOError: If the output file cannot be opened or written.
    """

    correct_hydrophobicity = []
    residues = ["I", "V", "L", "F", "M", "A", "W", "T"]

    for protein in proteinlist:
        percentages = ProteinAnalysis(str(protein.seq)).get_amino_acids_percent()

        hydro_percent = sum(percentages[x] for x in residues)

        if hydro_percent > minvalue and hydro_percent < maxvalue:
            correct_hydrophobicity.append(protein)

    SeqIO.write(correct_hydrophobicity, outputfile, "fasta")

    return(correct_hydrophobicity)

#filter for acidic residue percentage range
def acidicfilter(proteinlist, minvalue, maxvalue, outputfile):
    """
    Filters a list of SeqRecord objects based on the acidity of their protein sequence.
    Hydrophobicity is defined as the percentage of acidic residues [EDC]

    Args:
        proteinlist (list): A list of SeqRecord objects.
        minvalue (float): The minimum acidity of the protein sequence to include.
        maxvalue (float): The maximum acidity of the protein sequence to include.
        outputfile (str): The path to the output FASTA file.

    Returns:
        list: A list of filtered SeqRecord objects.

    Raises:
        IOError: If the output file cannot be opened or written.
    """


    correct_acidic = []
    residues = ["E", "D", "C"]

    for protein in proteinlist:
        percentages = ProteinAnalysis(str(protein.seq)).get_amino_acids_percent()

        acidic_percent = sum(percentages[x] for x in residues)

        if acidic_percent > minvalue and acidic_percent < maxvalue:
            correct_acidic.append(protein)

    SeqIO.write(correct_acidic, outputfile, "fasta")

    return(correct_acidic)
