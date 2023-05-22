from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sequence_processing import remove_illegal_chars
from sequence_processing import dictify

residue1 = "E"
residue2 = ["A", "S", "V", "F", "H", "T", "N", "I"]
residue3 = ["L", "V", "I", "F", "M", "T"]
lenientacidic = ["E", "D", "C"]
lenienthydrophobic = ["I", "V", "L", "F", "M", "W", "T", "A"]
acidic = ["E", "D"]
hydrophobic = ["I", "V", "L", "F", "M"]

def get_pi_values(proteinlist):
    """
    Calculate the isoelectric point (pI) values for a list of proteins.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of calculated pI values for the proteins.

    """

    pi_values = []

    for protein in proteinlist:
        pi_values.append(ProteinAnalysis(str(protein.seq)).isoelectric_point())

    return pi_values

def get_lengths(proteinlist):
    """
    Retrieve the lengths of proteins in a list.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of lengths (number of amino acids) for the proteins.

    """

    lengths = []

    for protein in proteinlist:
        lengths.append(len(str(protein.seq)))

    return lengths

def get_mws(proteinlist):
    """
    Calculate the molecular weights (MWs) of proteins in a list.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of calculated MWs (in kilodaltons) for the proteins.

    """

    weights = []

    for protein in proteinlist:
        weights.append(ProteinAnalysis(remove_illegal_chars(str(protein.seq))).molecular_weight() / 1000)

    return weights

def get_hydrophobicities(proteinlist):
    """
    Calculate the hydrophobicity scores for proteins in a list.

    The hydrophobicity score is calculated as the sum of percentages of hydrophobic amino acids (defined in the variable 'hydrophobic') in the protein sequence.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of calculated hydrophobicity scores for the proteins.

    """

    hydrophobicities = []

    for protein in proteinlist:
        percentages = ProteinAnalysis(str(protein.seq)).get_amino_acids_percent()

        hydrophobicities.append(sum(percentages[x] for x in hydrophobic) * 100)
    
    return hydrophobicities

def get_lenient_hydrophobicities(proteinlist):

    """
    Calculate the lenient hydrophobicity scores for proteins in a list.

    The lenient hydrophobicity score is calculated as the sum of percentages of lenient hydrophobic amino acids (defined in the variable 'lenienthydrophobic') in the protein sequence.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of calculated lenient hydrophobicity scores for the proteins.

    """

    hydrophobicities = []

    for protein in proteinlist:
        percentages = ProteinAnalysis(str(protein.seq)).get_amino_acids_percent()

        hydrophobicities.append(sum(percentages[x] for x in lenienthydrophobic) * 100)
    
    return hydrophobicities

def get_acidities(proteinlist):
    """
    Calculate the acidity scores for proteins in a list.

    The acidity score is calculated as the sum of percentages of acidic amino acids (defined in the variable 'acidic') in the protein sequence.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of calculated acidity scores for the proteins.

    """

    acidities = []

    for protein in proteinlist:
        percentages = ProteinAnalysis(str(protein.seq)).get_amino_acids_percent()

        acidities.append(sum(percentages[x] for x in acidic) * 100)
    
    return acidities

def get_lenient_acidities(proteinlist):
    """
    Calculate the lenient acidity scores for proteins in a list.

    The lenient acidity score is calculated as the sum of percentages of lenient acidic amino acids (defined in the variable 'lenientacidic') in the protein sequence.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of calculated lenient acidity scores for the proteins.

    """

    acidities = []

    for protein in proteinlist:
        percentages = ProteinAnalysis(str(protein.seq)).get_amino_acids_percent()

        acidities.append(sum(percentages[x] for x in lenientacidic) * 100)
    
    return acidities

def get_glypros(proteinlist):
    """
    Calculate the counts of glycine and proline residues for proteins in a list.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of counts of glycine and proline residues for the proteins.

    """

    glypros = []

    for protein in proteinlist:
        glypros.append(ProteinAnalysis(str(protein.seq)).count_amino_acids()["P"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["G"])
    
    return glypros

def get_ratios(proteinlist):
    """
    Calculate the acidic to basic amino acid ratios for proteins in a list.

    The ratio is calculated as the count of acidic amino acids (E and D) divided by the count of basic amino acids (K and R).

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of calculated acidic to basic amino acid ratios for the proteins.

    """

    ratios = []

    for protein in proteinlist:
        acidcount = ProteinAnalysis(str(protein.seq)).count_amino_acids()["E"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["D"]
        basecount = ProteinAnalysis(str(protein.seq)).count_amino_acids()["K"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["R"]

        try:
            ratios.append(acidcount / basecount)
        except ZeroDivisionError:
            ratios.append(0)
    
    return ratios

def get_lenient_ratios(proteinlist):

    """
    Calculate the lenient acidic to basic amino acid ratios for proteins in a list.

    The ratio is calculated as the count of lenient acidic amino acids (E, D, and C) divided by the count of lenient basic amino acids (K, R, and H).

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of calculated lenient acidic to basic amino acid ratios for the proteins.

    """

    ratios = []

    for protein in proteinlist:
        acidcount = ProteinAnalysis(str(protein.seq)).count_amino_acids()["E"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["D"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["C"]
        basecount = ProteinAnalysis(str(protein.seq)).count_amino_acids()["K"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["R"] + ProteinAnalysis(str(protein.seq)).count_amino_acids()["H"]

        try:
            ratios.append(acidcount / basecount)
        except ZeroDivisionError:
            ratios.append(0)
    
    return ratios

def get_esi_glypros(proteinlist):

    """
    Calculate the counts of glycine and proline residues within the ESI motif for proteins in a list.

    The ESI motif is defined as the presence of specific residues in a specific pattern.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of counts of glycine and proline residues within the ESI motif for the proteins.

    """

    esi_glypros = []

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
                            esi_glypros.append( ProteinAnalysis(strseq[position:position + 20]).count_amino_acids()["P"] + ProteinAnalysis(strseq[position:position + 20]).count_amino_acids()["G"])
                            break
            except KeyError:
                break
    
    return esi_glypros

def get_esi_distances(proteinlist):

    """
    Calculate the distances of the ESI motif from the C-terminus for proteins in a list.

    The ESI motif is defined as the presence of specific residues in a specific pattern.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of calculated distances of the ESI motif from the C-terminus for the proteins.

    """

    esi_distances = []

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+1] in residue2:
                        if sequence[i+2] in residue3:
                            esi_distances.append(len(str(protein.seq)) - i + 1)
                            break
            except KeyError:
                break
    
    return esi_distances

def get_esi_aspglus(proteinlist):
    """
    Calculate the counts of acidic residues within a specific region around the ESI motif for proteins in a list.

    The ESI motif is defined as the presence of specific residues in a specific pattern.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of counts of acidic residues within the region around the ESI motif for the proteins.

    """

    esi_aspglus = []

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
                                    esi_aspglus.append(ProteinAnalysis(strseq[position-6:position + 22]).count_amino_acids()["E"] + ProteinAnalysis(strseq[position:position-6 + 22]).count_amino_acids()["D"])
                                    break
                                except IndexError:
                                    break  
            except KeyError:
                break
        
        return esi_aspglus

def get_esi_acidic(proteinlist):
    """
    Determine if the ESI motif in proteins from a list contains acidic residues.

    The ESI motif is defined as the presence of specific residues in a specific pattern.

    Args:
        proteinlist (list): List of SeqRecord objects representing the proteins.

    Returns:
        list: List of binary values indicating the presence (1) or absence (0) of acidic residues within the ESI motif for the proteins.

    """

    esi_acidic = []

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
                                        esi_acidic.append(1)
                                        break
                                    else:
                                        esi_acidic.append(0)
                                        break
                            except IndexError:
                                break
            except KeyError:
                break
    
    return esi_acidic