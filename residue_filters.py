from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sequence_processing import remove_duplicates
import re

#filter for glycine + proline count range
def glyprofilter(proteinlist, mincount, maxcount, outputfile):

    correct_glypro = []
    
    for protein in proteinlist:
        count = ProteinAnalysis(str(protein.seq)).count_amino_acids()["P"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["G"]

        if count >= mincount and count <= maxcount:
            correct_glypro.append(protein)

    SeqIO.write(correct_glypro, outputfile, "fasta")

    return(correct_glypro)

#filter for min acidic:basic ratio
def chargeratiofilter(proteinlist, minratio, outputfile):

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

    correct_motifs = []

    for protein in proteinlist:

        z_sequence = re.sub("(EE|DD|DE|ED)", "Z", str(protein.seq))

        ed_list = re.findall("Z", z_sequence)
        
        if len(ed_list) >= 3:
            correct_motifs.append(protein)

    SeqIO.write(correct_motifs, outputfile, "fasta")

    return(correct_motifs)