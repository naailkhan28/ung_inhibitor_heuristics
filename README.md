# Ung Inhibitor Heuristics

This GitHub repo contains data and scripts from our paper *A multimodal approach towards genomic identification of protein inhibitors of uracil-DNA glycosylase*.

## Scripts

The file `filter_runs.py` contains an example heuristics filter run. The broad workflow for a heuristics run is as follows:

1. Read in data from a single list of proteins (`single_inputs.py`) or a whole genome or list of genomes (`genome_inputs.py`)

2. Remove duplicates and trim translated sequences to the first start codon (`sequence_processing.py`)

3. Filter sequences by acidity and hydrophobicity (`general_filters.py`)

4. Filter based on glycine/proline residues or ratio of acidic and basic residues (`residue_filters.py`)

5. Filter based on presence of ESI motif and surrounding residues (either `esi_filters.py` or `lenient_esi_filters.py` which follow the lenient definition of the ESI motif).

### Additional Scripts

`plot_histograms.py` contains scripts for plotting the distribution of lengths, molecular weights, isoelectric points, and various other physical properties across a set of proteins. This makes use of `sequence_analysis.py` which contains functions to get these properties for a single protein sequence.

`plot_hit_rates.py` is a simple dummy plot script to plot a bar chart from a list of floating point values

`p56_filters.py` contains functions for filtering protein sequences based on p56-type Ung inhibitor motifs (EXXYG and FXDSY). These are similar to the ESI-motif filters for Ugi and SAUGI-type UngIns but were not presented or used in the manuscript.

`window_analysis.py` and `window_filters.py` contain functions to filter proteins based on sliding-window metrics such as molecular weight and pI of discrete stretches of residues. These were explored during heuristics development to measure characteristics of putative β-strands containing ESI motifs, for example. This was not presented in the manuscript but scripts are provided here for reference.

## Data

There are six data files:

 1. `All_Ugis+SaUgis.fasta` contains protein sequences for all Ugi and SaUgi-type Ung inhibitors, used to make curated MSAs
 2. `supplementary_dna_sequences.fasta` contains DNA sequences from Supplementary Section 1 of the Supplementary Information. 
 3. `supplementary_protein_sequences.fasta` contains protein sequences from Supplementary Section 1 of the Supplementary Information. 
 4. `supplementary_table_s6_heuristics_matches.fasta` contains protein sequences from Supplementary Table S6 of the Supplementary Information, corresponding to selected heuristics hit sequences. 
 5. `Filter Tests.xlsx` contains information and data from testing and tuning of heuristics filters, as well as information about each filter script.
 6. `s6_genome.fasta` is an example genome from *Staphylococcus* phage S6. 