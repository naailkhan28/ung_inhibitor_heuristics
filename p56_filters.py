from Bio import SeqIO
from sequence_processing import dictify

#filter for p56 EXXYG-like motifs
def lenient_exfinder(proteinlist, outputfile):

    correct_motif = []

    residue1 = ["E", "G", "S", "H", "C"]
    residue3 = ["L", "M", "V", "I"]
    residue4 = ["Y", "V", "I"]
    residue5 = "G"

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] in residue1:
                    if sequence[i+2] in residue3:
                        if sequence[i+3] in residue4:
                            if sequence[i+4] == residue5:        
                                correct_motif.append(protein)
                                break
            except KeyError:
                break

    SeqIO.write(correct_motif, outputfile, "fasta")

    return(correct_motif)

#filter for FXDSY-like motifs
def lenient_fxfinder(proteinlist, outputfile):

    correct_motif = []

    residue1 = ["F", "K", "L", "E", "D", "N", "G", "Y"]
    residue3 = ["D", "E", "Q"]
    residue4 = ["S", "T", "L", "F", "W"]
    residue5 = ["Y", "W", "C"]

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] in residue1:
                    if sequence[i+2] in residue3:
                        if sequence[i+3] in residue4:
                            if sequence[i+4] in residue5:        
                                correct_motif.append(protein)
                                break
            except KeyError:
                break

    SeqIO.write(correct_motif, outputfile, "fasta")

    return(correct_motif)

def lenient_motifdistance(proteinlist, mindistance, maxdistance, outputfile):

    correct_distance = []

    fx_residue1 = ["F", "K", "L", "E", "D", "N", "G", "Y"]
    fx_residue3 = ["D", "E", "Q"]
    fx_residue4 = ["S", "T", "L", "F", "W"]
    fx_residue5 = ["Y", "W", "C"]

    ex_residue1 = ["E", "G", "S", "H", "C"]
    ex_residue3 = ["L", "M", "V", "I"]
    ex_residue4 = ["Y", "V", "I"]
    ex_residue5 = "G"
    
    for protein in proteinlist:
        
        sequence = dictify(str(protein.seq))

        fx = 0
        ex = 0

        for i in sequence:

            try:
                if sequence[i] in fx_residue1:
                    if sequence[i+2] in fx_residue3:
                        if sequence[i+3] in fx_residue4:
                            if sequence[i+4] in fx_residue5:        
                                fx = i
                if sequence[i] in ex_residue1:
                    if sequence[i+2] in ex_residue3:
                        if sequence[i+3] in ex_residue4:
                            if sequence[i+4] == ex_residue5:        
                                ex = i
                                if ex != 0 and fx != 0:
                                    if ex - fx >= mindistance and ex - fx <= maxdistance:
                                        #print(ex, fx, fx-ex)
                                        correct_distance.append(protein)
                                        break
            
            except KeyError:
                break
        
    SeqIO.write(correct_distance, outputfile, "fasta")

    return(correct_distance)

    #filter for p56 EXXYG-like motifs
def strict_exfinder(proteinlist, outputfile):

    correct_motif = []

    residue1 = ["E", "G", "S"]
    residue3 = ["L", "M", "V", "I"]
    residue4 = ["Y", "V", "I"]
    residue5 = "G"

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] in residue1:
                    if sequence[i+2] in residue3:
                        if sequence[i+3] in residue4:
                            if sequence[i+4] == residue5:        
                                correct_motif.append(protein)
                                break
            except KeyError:
                break

    SeqIO.write(correct_motif, outputfile, "fasta")

    return(correct_motif)

#filter for FXDSY-like motifs
def strict_fxfinder(proteinlist, outputfile):

    correct_motif = []

    residue1 = "F"
    residue3 = "D"
    residue4 = "S"
    residue5 = "Y"

    for protein in proteinlist:

        sequence = dictify(str(protein.seq))

        for i in sequence:
            try:
                if sequence[i] == residue1:
                    if sequence[i+2] == residue3:
                        if sequence[i+3] == residue4:
                            if sequence[i+4] == residue5:        
                                correct_motif.append(protein)
                                break
            except KeyError:
                break

    SeqIO.write(correct_motif, outputfile, "fasta")

    return(correct_motif)

def strict_motifdistance(proteinlist, mindistance, maxdistance, outputfile):

    correct_distance = []

    fx_residue1 = "F"
    fx_residue3 = "D"
    fx_residue4 = "S"
    fx_residue5 = "Y"

    ex_residue1 = ["E", "G", "S"]
    ex_residue3 = ["L", "M", "V", "I"]
    ex_residue4 = ["Y", "V", "I"]
    ex_residue5 = "G"
    
    for protein in proteinlist:
        
        sequence = dictify(str(protein.seq))

        fx = 0
        ex = 0

        for i in sequence:

            try:
                if sequence[i] == fx_residue1:
                    if sequence[i+2] == fx_residue3:
                        if sequence[i+3] == fx_residue4:
                            if sequence[i+4] == fx_residue5:        
                                fx = i
                if sequence[i] in ex_residue1:
                    if sequence[i+2] in ex_residue3:
                        if sequence[i+3] in ex_residue4:
                            if sequence[i+4] == ex_residue5:        
                                ex = i
                                if ex != 0 and fx != 0:
                                    if ex - fx >= mindistance and ex - fx <= maxdistance:
                                        #print(ex, fx, fx-ex)
                                        correct_distance.append(protein)
                                        break
            
            except KeyError:
                break
        
    SeqIO.write(correct_distance, outputfile, "fasta")

    return(correct_distance)