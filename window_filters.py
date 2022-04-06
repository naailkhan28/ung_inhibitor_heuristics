from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from sequence_processing import remove_illegal_chars
from sequence_processing import dictify

residue1 = "E"
residue2 = ["A", "S", "V", "F", "H", "T", "N", "I"]
residue3 = ["L", "V", "I", "F", "M", "T"]
lenienthydrophobic = ["I", "V", "L", "F", "M", "W", "T", "A"]

def esi_strand_mw(proteinlist, minmw, maxmw, outputfile):

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
                                if ProteinAnalysis(remove_illegal_chars(strand)).molecular_weight() >= minmw and ProteinAnalysis(remove_illegal_chars(strand)).molecular_weight() <= maxmw:
                                    esi_mws.append(protein)
                                break
            except KeyError:
                break
    SeqIO.write(esi_mws, outputfile, "fasta")
    return esi_mws

def esi_strand_pi(proteinlist, minpi, maxpi, outputfile):

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
                                if ProteinAnalysis(remove_illegal_chars(strand)).isoelectric_point() >= minpi and ProteinAnalysis(remove_illegal_chars(strand)).isoelectric_point() <= maxpi:
                                    esi_pis.append(protein)
                                break
            except KeyError:
                break
    SeqIO.write(esi_pis, outputfile, "fasta")
    return esi_pis

def esi_strand_hydrophobicity(proteinlist, minhydro, maxhydro, outputfile):

    esi_hydros = []

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

                                strand_percent = sum(percentages[x] for x in lenienthydrophobic) * 100
                                if strand_percent >= minhydro and strand_percent <= maxhydro:
                                    esi_hydros.append(protein)
                                break
            except KeyError:
                break
    SeqIO.write(esi_hydros, outputfile, "fasta")
    return esi_hydros

def lenient_esi_strand_mw(proteinlist, minmw, maxmw, outputfile):

    esi_mws = []

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+2] in residue3:
                        if i > 1 and i + 4 < len(str(protein.seq)):
                            strand = str(protein.seq)[i-2:i+5]
                            if ProteinAnalysis(remove_illegal_chars(strand)).molecular_weight() >= minmw and ProteinAnalysis(remove_illegal_chars(strand)).molecular_weight() <= maxmw:
                                esi_mws.append(protein)
                            break
            except KeyError:
                break
    SeqIO.write(esi_mws, outputfile, "fasta")
    return esi_mws

def lenient_esi_strand_pi(proteinlist, minpi, maxpi, outputfile):

    esi_pis = []

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+2] in residue3:
                        if i > 1 and i + 4 < len(str(protein.seq)):
                            strand = str(protein.seq)[i-2:i+5]
                            if ProteinAnalysis(remove_illegal_chars(strand)).isoelectric_point() >= minpi and ProteinAnalysis(remove_illegal_chars(strand)).isoelectric_point() <= maxpi:
                                esi_pis.append(protein)
                            break
            except KeyError:
                break
    SeqIO.write(esi_pis, outputfile, "fasta")
    return esi_pis

def lenient_esi_strand_hydrophobicity(proteinlist, minhydro, maxhydro, outputfile):

    esi_hydros = []

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+2] in residue3:
                        if i > 1 and i + 4 < len(str(protein.seq)):
                            strand = str(protein.seq)[i-2:i+5]
                            percentages = ProteinAnalysis(strand).get_amino_acids_percent()

                            strand_percent = sum(percentages[x] for x in lenienthydrophobic) * 100
                            if strand_percent >= minhydro and strand_percent <= maxhydro:
                                esi_hydros.append(protein)
                            break
            except KeyError:
                break
    SeqIO.write(esi_hydros, outputfile, "fasta")
    return esi_hydros