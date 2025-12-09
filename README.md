# Analysis-Tool-for-Local-Association-of-Neuronal-Transcript-Expression-ATLANTE-
ATLANTE identifies which brain regions show unusually high expression for a custome list of genes of interest (human). Uses Monte Carlo simulation to compare your gene list against 193 human brain regions from the Human Protein Atlas.

Overview
ATLANTE is a Python-based analysis pipeline designed to identify brain regions with statistically significant enrichment of highly expressed genes from user-defined gene lists. The tool uses Monte Carlo simulation to compare observed gene expression patterns against reference distributions, enabling researchers to map genetic findings onto discrete neuroanatomical targets.

Key Features
- Monte Carlo simulation (default: 1,000,000 iterations) with Bonferroni correction
- Comprehensive coverage: Analyzes 193 discrete human brain regions
- Utilizes ~20,162 genes from the Human Protein Atlas
- Successfully applied to Major Depression and Autism Spectrum Disorder gene sets
- Works with any user-defined gene list

Installation
Requirements
Python 3.8+

Dependencies
pip install pandas numpy scipy matplotlib plotly networkx

Specific versions used in development
pandas 2.2.2
numpy 2.3.4
scipy 1.13.1
matplotlib 3.9.2
plotly 6.3.1
networkx 3.5

Data Dependencies
Required Data File:
human_protein_atlas_expression.csv

This file contains RNA transcript counts for 193 brain regions from the Human Protein Atlas (version 24.0 or later).

Format: 
Rows=Individual genes (20,162+ entries)
Columns=Brain regions (193 regions)
Values=Normalized transcript expression levels
Source=Download from Human Protein Atlas (www.proteinatlas.org)
Note - Due to size and licensing, this file is not included in the repository. Users must download it independently.

Citation
If you use ATLANTE in your research, please cite:
Odierna, G.L., Sharpley, C.F., Bitsika, V., Evans, I.D., & Vessey, K.A. (2025). 
Novel Gene-Informed Regional Brain Targets for Clinical Screening for Major Depression. 
Neurology International, 17(6). doi:10.3390/neurolint17060096
