from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sequence_processing import remove_illegal_chars
from sequence_processing import dictify

residue1 = "E"
residue2 = ["A", "S", "V", "F", "H", "T", "N", "I"]
residue3 = ["L", "V", "I", "F", "M", "T"]
lenienthydrophobic = ["I", "V", "L", "F", "M", "W", "T", "A"]

def get_esi_mws(proteinlist):

    esi_mws = []

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+1] in residue2:
                        if sequence[i+2] in residue3:
                            if i > 1 and i + 4 < len(str(protein.seq)):
                                strand = str(protein.seq)[i-2:i+5]
                                esi_mws.append(ProteinAnalysis(remove_illegal_chars(strand)).molecular_weight())
                                break
            except KeyError:
                break
    
    return esi_mws

def get_esi_pi(proteinlist):

    esi_pis = []

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+1] in residue2:
                        if sequence[i+2] in residue3:
                            if i > 1 and i + 4 < len(str(protein.seq)):
                                strand = str(protein.seq)[i-2:i+5]
                                esi_pis.append(ProteinAnalysis(remove_illegal_chars(strand)).isoelectric_point())
                                break
            except KeyError:
                break
    
    return esi_pis

def get_eight_residue_pi(proteinlist):

    pis = []

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if i > 1 and i + 4 < len(str(protein.seq)):
                    strand = str(protein.seq)[i-2:i+5]
                    pis.append(ProteinAnalysis(remove_illegal_chars(strand)).isoelectric_point())
            except KeyError:
                break
    
    return pis

def get_eight_residue_mw(proteinlist):

    residues_mws = []

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if i > 1 and i + 4 < len(str(protein.seq)):
                    strand = str(protein.seq)[i-2:i+5]
                    residues_mws.append(ProteinAnalysis(remove_illegal_chars(strand)).molecular_weight())
            except KeyError:
                break
    
    return residues_mws

def get_esi_hydrophobicity(proteinlist):

    hydros = []

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+1] in residue2:
                        if sequence[i+2] in residue3:
                            if i > 1 and i + 4 < len(str(protein.seq)):
                                strand = str(protein.seq)[i-2:i+5]
                                percentages = ProteinAnalysis(strand).get_amino_acids_percent()

                                hydros.append(sum(percentages[x] for x in lenienthydrophobic) * 100)
                                break
            except KeyError:
                break
    
    return hydros

def get_eight_residue_hydrophobicity(proteinlist):

    hydros = []

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if i > 1 and i + 4 < len(str(protein.seq)):
                    strand = str(protein.seq)[i-2:i+5]
                    percentages = ProteinAnalysis(strand).get_amino_acids_percent()

                    hydros.append(sum(percentages[x] for x in lenienthydrophobic) * 100)
            except KeyError:
                break
    
    return hydros