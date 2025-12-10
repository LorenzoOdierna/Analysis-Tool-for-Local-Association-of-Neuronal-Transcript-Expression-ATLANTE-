# Analysis-Tool-for-Local-Association-of-Neuronal-Transcript-Expression-ATLANTE-
ATLANTE identifies which brain regions show unusually high expression for a custome list of genes of interest (human). Uses Monte Carlo simulation to compare your gene list against 193 human brain regions from the Human Protein Atlas.


**Installation**

Requirements

Python 3.8+


**Dependencies**

pandas numpy  



**Data Dependencies**

Required Data File:

rna_brain_hpa.tsv
Source=Download from Human Protein Atlas (www.proteinatlas.org)

This file contains RNA transcript counts for 193 brain regions from the Human Protein Atlas (version 24.0 or later).

Format: 

Rows=Individual genes (20,162+ entries)

Columns=Brain regions (193 regions)

Values=Normalized transcript expression levels

Note - Due to size and licensing, this file is not included in the repository. Users must download it independently.


**Citation**

If you use ATLANTE in your research, please cite:

Odierna, G.L., Sharpley, C.F., Bitsika, V., Evans, I.D., & Vessey, K.A. (2025). 
Novel Gene-Informed Regional Brain Targets for Clinical Screening for Major Depression. 
Neurology International, 17(6). doi:10.3390/neurolint17060096
