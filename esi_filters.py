from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sequence_processing import dictify

#filter for ESI-like motifs
def esifinder(proteinlist, outputfile):

    correct_motif = []

    residue1 = "E"
    residue2 = ["A", "S", "V", "F", "H", "T", "N", "I"]
    residue3 = ["L", "V", "I", "F", "M", "T"]

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+1] in residue2:
                        if sequence[i+2] in residue3:
                            correct_motif.append(protein)
                            break
            except KeyError:
                break

    SeqIO.write(correct_motif, outputfile, "fasta")

    return(correct_motif)

#filter for length between ESI-like motif and stop codon
def esidistance(proteinlist, mindistance, maxdistance, outputfile):

    correct_distance= []
    
    residue1 = "E"
    residue2 = ["A", "S", "V", "F", "H", "T", "N", "I"]
    residue3 = ["L", "V", "I", "F", "M", "T"]

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+1] in residue2:
                        if sequence[i+2] in residue3:
                            distance = len(str(protein.seq)) - (i + 1)
                            if distance > mindistance and distance < maxdistance:
                                correct_distance.append(protein)
                                break
            except KeyError:
                break

    SeqIO.write(correct_distance, outputfile, "fasta")

    return(correct_distance)

#filter for glycine and proline residues around ESI-like motif
def esiresidues_glypro(proteinlist, outputfile):

    correct_esiresidues_glypro = []
    
    residue1 = "E"
    residue2 = ["A", "S", "V", "F", "H", "T", "N", "I"]
    residue3 = ["L", "V", "I", "F", "M", "T"]

    for protein in proteinlist:

        position = -1

        strseq = str(protein.seq)

        sequence = dictify(strseq)

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+1] in residue2:
                        if sequence[i+2] in residue3:
                            position = i
                            try:
                                gpcount = ProteinAnalysis(strseq[position:position + 20]).count_amino_acids()["P"] + ProteinAnalysis(strseq[position:position + 20]).count_amino_acids()["G"]
                                if gpcount == 2 or gpcount == 3:
                                    correct_esiresidues_glypro.append(protein)
                                    break
                            except IndexError:
                                    break
            except KeyError:
                break
	
    SeqIO.write(correct_esiresidues_glypro, outputfile, "fasta")

    return(correct_esiresidues_glypro)

#filter for acidic residues around ESI-like motif
def esiresidues_acidic(proteinlist, outputfile):


    correct_esiresidues_acidic = []
    
    residue1 = "E"
    residue2 = ["A", "S", "V", "F", "H", "T", "N", "I"]
    residue3 = ["L", "V", "I", "F", "M", "T"]
    acidic = ["E", "D"]

    for protein in proteinlist:

        position = -1

        strseq = str(protein.seq)

        sequence = dictify(strseq)

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+1] in residue2:
                        if sequence[i+2] in residue3:
                            position = i
                            try:
                                if position != -1:
                                    if strseq[position+7] in acidic or strseq[position+8] in acidic:
                                        correct_esiresidues_acidic.append(protein)
                                        break
                            except IndexError:
                                break
            except KeyError:
                break

    SeqIO.write(correct_esiresidues_acidic, outputfile, "fasta")

    return(correct_esiresidues_acidic)

#filter for aspartate and glutamate residues around ESI-like motif
def esiresidues_aspglu(proteinlist, outputfile):



    correct_esiresidues_aspglu = []
    
    residue1 = "E"
    residue2 = ["A", "S", "V", "F", "H", "T", "N", "I"]
    residue3 = ["L", "V", "I", "F", "M", "T"]

    for protein in proteinlist:

        position = -1

        strseq = str(protein.seq)

        sequence = dictify(strseq)

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+1] in residue2:
                        if sequence[i+2] in residue3:
                            position = i
                            if position != -1:
                                try:
                                    edcount = ProteinAnalysis(strseq[position-6:position + 22]).count_amino_acids()["E"] + ProteinAnalysis(strseq[position:position-6 + 22]).count_amino_acids()["D"]
                                    if edcount >= 7:
                                        correct_esiresidues_aspglu.append(protein)
                                        break
                                except IndexError:
                                    break  
            except KeyError:
                break

    SeqIO.write(correct_esiresidues_aspglu, outputfile, "fasta")

    return(correct_esiresidues_aspglu)

def allesifilters(proteinlist, mindistance, maxdistance, outputfile):
    output= []
    
    residue1 = "E"
    residue2 = ["A", "S", "V", "F", "H", "T", "N", "I"]
    residue3 = ["L", "V", "I", "F", "M", "T"]
    acidic = ["E", "D"]

    for protein in proteinlist:

        position = -1

        strseq = str(protein.seq)

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+1] in residue2:
                        if sequence[i+2] in residue3:
                            distance = len(str(protein.seq)) - (i + 1)
                            if distance > mindistance and distance < maxdistance:
                                position = i
                                try:
                                    gpcount = ProteinAnalysis(strseq[position:position + 20]).count_amino_acids()["P"] + ProteinAnalysis(strseq[position:position + 20]).count_amino_acids()["G"]
                                    if gpcount == 2 or gpcount == 3:
                                        if strseq[position+7] in acidic or strseq[position+8] in acidic:
                                            edcount = ProteinAnalysis(strseq[position-6:position + 22]).count_amino_acids()["E"] + ProteinAnalysis(strseq[position:position-6 + 22]).count_amino_acids()["D"]
                                            if edcount >= 7:
                                                output.append(protein)
                                                break
                                except IndexError:
                                    break
            except KeyError:
                break

    SeqIO.write(output, outputfile, "fasta")

    return(output)