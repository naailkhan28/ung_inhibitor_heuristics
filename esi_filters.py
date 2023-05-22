from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sequence_processing import dictify

#filter for ESI-like motifs
def esifinder(proteinlist, outputfile):
    """
    Filter a list of proteins for ESI-like motifs.

    The function searches for a specific motif pattern in the proteins where:
    - The first residue is "E"
    - The second residue can be any of ["A", "S", "V", "F", "H", "T", "N", "I"]
    - The third residue can be any of ["L", "V", "I", "F", "M", "T"]

    Proteins containing the ESI-like motif will be added to the 'correct_motif' list.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.
        outputfile (str): Path to the output file where proteins with the ESI-like motif will be written in FASTA format.

    Returns:
        list: List of SeqRecord objects representing the proteins containing the ESI-like motif.

    """

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
    """
    Filter a list of proteins based on the distance between an ESI-like motif and the stop codon.

    The function searches for proteins that contain an ESI-like motif pattern, where:
    - The first residue is "E"
    - The second residue can be any of ["A", "S", "V", "F", "H", "T", "N", "I"]
    - The third residue can be any of ["L", "V", "I", "F", "M", "T"]

    Proteins are filtered based on the distance between the ESI-like motif and the stop codon,
    where the distance is greater than 'mindistance' and less than 'maxdistance'.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.
        mindistance (int): Minimum distance between the ESI-like motif and the stop codon.
        maxdistance (int): Maximum distance between the ESI-like motif and the stop codon.
        outputfile (str): Path to the output file where filtered proteins will be written in FASTA format.

    Returns:
        list: List of SeqRecord objects representing the proteins that meet the distance criteria.

    """
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
    """
    Filter a list of proteins based on the presence of glycine and proline residues around an ESI-like motif.

    The function searches for proteins that contain an ESI-like motif pattern, where:
    - The first residue is "E"
    - The second residue can be any of ["A", "S", "V", "F", "H", "T", "N", "I"]
    - The third residue can be any of ["L", "V", "I", "F", "M", "T"]

    Proteins are filtered based on the presence of glycine and proline residues within a 20 amino acid window around the ESI-like motif.
    If there are exactly 2 or 3 glycine and proline residues combined, the protein is considered to meet the criteria.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.
        outputfile (str): Path to the output file where filtered proteins will be written in FASTA format.

    Returns:
        list: List of SeqRecord objects representing the proteins that meet the glycine and proline residues criteria.

    """

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
    """
    Filter a list of proteins based on the presence of acidic residues around an ESI-like motif.

    The function searches for proteins that contain an ESI-like motif pattern, where:
    - The first residue is "E"
    - The second residue can be any of ["A", "S", "V", "F", "H", "T", "N", "I"]
    - The third residue can be any of ["L", "V", "I", "F", "M", "T"]

    Proteins are filtered based on the presence of acidic residues (E or D) at positions +7 or +8 from the ESI-like motif.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.
        outputfile (str): Path to the output file where filtered proteins will be written in FASTA format.

    Returns:
        list: List of SeqRecord objects representing the proteins that meet the acidic residues criteria.

    """


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

    """
    Filter a list of proteins based on the presence of aspartate and glutamate residues around an ESI-like motif.

    The function searches for proteins that contain an ESI-like motif pattern, where:
    - The first residue is "E"
    - The second residue can be any of ["A", "S", "V", "F", "H", "T", "N", "I"]
    - The third residue can be any of ["L", "V", "I", "F", "M", "T"]

    Proteins are filtered based on the presence of aspartate (D) and glutamate (E) residues within a 28 amino acid window around the ESI-like motif.
    If there are 7 or more aspartate and glutamate residues combined, the protein is considered to meet the criteria.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.
        outputfile (str): Path to the output file where filtered proteins will be written in FASTA format.

    Returns:
        list: List of SeqRecord objects representing the proteins that meet the aspartate and glutamate residues criteria.

    """

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
    """
    Filter a list of proteins based on multiple criteria related to an ESI-like motif.

    The function filters proteins based on the following criteria:
    - Presence of an ESI-like motif pattern, where:
        - The first residue is "E"
        - The second residue can be any of ["A", "S", "V", "F", "H", "T", "N", "I"]
        - The third residue can be any of ["L", "V", "I", "F", "M", "T"]
    - Distance between the ESI-like motif and the stop codon of the protein, where:
        - The distance is greater than `mindistance` and less than `maxdistance`
    - Presence of glycine (G) and proline (P) residues within a 20 amino acid window around the ESI-like motif, where:
        - The count of glycine and proline residues is either 2 or 3
    - Presence of acidic residues (E or D) at positions +7 or +8 from the ESI-like motif
    - Presence of aspartate (D) and glutamate (E) residues within a 28 amino acid window around the ESI-like motif, where:
        - The count of aspartate and glutamate residues combined is 7 or more

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.
        mindistance (int): Minimum distance between the ESI-like motif and the stop codon.
        maxdistance (int): Maximum distance between the ESI-like motif and the stop codon.
        outputfile (str): Path to the output file where filtered proteins will be written in FASTA format.

    Returns:
        list: List of SeqRecord objects representing the proteins that meet all the filtering criteria.

    """
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