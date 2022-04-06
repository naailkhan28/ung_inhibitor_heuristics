import single_inputs
import genome_inputs
import sequence_processing
import general_filters
import residue_filters
import esi_filters
import lenient_esi_filters
import window_filters
import p56_filters

table = 11
min_length = 40
start = ["M", "V", "I", "L"]

proteins = genome_inputs.read_genomes_bulk("s6_genome.fasta", table, min_length, "1_proteins.fasta")

trimmed = sequence_processing.trim_residues(proteins, start, "2_trimmed.fasta")

hydrofilter = general_filters.hydrophobicfilter(trimmed, 0.36, 0.48, "3_hydrofilter.fasta")

acidfilter = general_filters.acidicfilter(hydrofilter, 0.12, 0.24, "4_acidfilter.fasta")

glypro = residue_filters.glyprofilter(acidfilter, 5, 10, "5_glypro.fasta")

ratio = residue_filters.chargeratiofilter(glypro, 1.0, "6_ratio.fasta")

esi = lenient_esi_filters.allesifilters(ratio, 50, 120, "7_esifilters.fasta")