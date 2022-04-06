from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#filter below max length and molecular weight
def lengthfilter(proteinlist, maxlength, maxmW, outputfile):

    correct_size = []

    for protein in proteinlist:
        if ProteinAnalysis(str(protein.seq)).molecular_weight() < maxmW and len(str(protein.seq)) < maxlength:
            correct_size.append(protein)

    SeqIO.write(correct_size, outputfile, "fasta")

    return(correct_size)

#filter for correct pI range
def pifilter(proteinlist, minpi, maxpi, outputfile):

    correct_pi = []

    for protein in proteinlist:
        if ProteinAnalysis(str(protein.seq)).isoelectric_point() <= maxpi and ProteinAnalysis(str(protein.seq)).isoelectric_point() >= minpi:
            correct_pi.append(protein)

    SeqIO.write(correct_pi, outputfile, "fasta")

    return(correct_pi)

#filter for hydrophobic residue percentage range
def hydrophobicfilter(proteinlist, minvalue, maxvalue, outputfile):

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

    correct_acidic = []
    residues = ["E", "D", "C"]

    for protein in proteinlist:
        percentages = ProteinAnalysis(str(protein.seq)).get_amino_acids_percent()

        acidic_percent = sum(percentages[x] for x in residues)

        if acidic_percent > minvalue and acidic_percent < maxvalue:
            correct_acidic.append(protein)

    SeqIO.write(correct_acidic, outputfile, "fasta")

    return(correct_acidic)
