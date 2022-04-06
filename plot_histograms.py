import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sequence_analysis
import window_analysis
from single_inputs import read_proteins
import numpy as np

protelomerases = read_proteins("uncharacterized_sequences.fasta")

lengths = sequence_analysis.get_lengths(protelomerases)

# lengths = sequence_analysis.get_lengths(p56)
# weights = sequence_analysis.get_mws(p56)
# hydrophobic = sequence_analysis.get_hydrophobicities(p56)
# acidity = sequence_analysis.get_acidities(p56)
# pi = sequence_analysis.get_pi_values(p56)
# glypro = sequence_analysis.get_glypros(p56)
# ratio = sequence_analysis.get_ratios(p56)

# lenient_acidity = sequence_analysis.get_lenient_acidities(p56)
# lenient_hydrophobic = sequence_analysis.get_lenient_hydrophobicities(p56)
# lenient_ratio = sequence_analysis.get_lenient_ratios(p56)

# esi_mws = window_analysis.get_esi_mws(p56)
# esi_pis = window_analysis.get_esi_pi(p56)
# esi_hydros = window_analysis.get_esi_hydrophobicity(p56)
# eight_mw = window_analysis.get_eight_residue_mw(p56)
# eight_pis = window_analysis.get_eight_residue_pi(p56)
# eight_hydros = window_analysis.get_eight_residue_hydrophobicity(p56)

fig = plt.hist(lengths)

# axs[0].hist(hydrophobic,bins=np.arange(min(hydrophobic), max(hydrophobic) + 2, 2), edgecolor="black", color="#FCF186")
# axs[0].set_title('Hydrophobic Residues [%]')

# axs[1].hist(acidity,bins=np.arange(min(acidity), max(acidity) + 2, 2), edgecolor="black", color="#74AFB7")
# axs[1].set_title('Acidic Residues [%]')

# for ax in axs.flat:
#     ax.set(ylabel='No. of Sequences')

plt.show()