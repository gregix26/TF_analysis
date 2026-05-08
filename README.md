# TF_analysis
## Choice 1 (safe): HOMER
Scan a defined window around SNP for human TF motifs with HOMER. 

### Setup
Dependencies: homer, r-base, samtools, r-essentials, bioconductor-deseq2, bioconductor-edger, bioconda-homer, bedtools

Easiest way to get HOMER is via conda. 

Within HOMER, create custom motif file if you want the analysis to be targeted.
First, create a txt file of TF families you are interested in. 
Compile major human TF in the motif directory: grep -Ei 'ctcf|ap1|jun|fos|nfkb|stat|fox|gata|sox|pax|ets|runx|smad|ceb|klf|e2f|irf|myc|max'     <(ls ~/.conda/envs/HOMER/share/homer-4.10-0/motifs/)     > matched_files.txt
Create a new motif file matching those you are interested in: while read f; do cat ~/.conda/envs/HOMER/share/homer-4.10-0/motifs/$f >> human_only.motifs; done < matched_files.txt 

### Run
In this version of HOMER, a 20bp window is extracted around the SNP from reference FASTA file for reference and alternate allele. HOMER scans the sequences for motifs. The R script calculates the delta value of REF - ALT motif score which informs about TF binding motif gained, lost, weaked or strenghtened. 

## Choice 2 (advanced)
Make use of CRESted and TF-MINDI for the prediciton of motif syntax disruption by variant in a cell-type-specific manner

### Setup

Pip install TF-MINDI and CRESted[motif] (with Tensorflow and cuda), modisco (motif discovery)
Dependendies: pybedtools pysam

There's two human brain pre-trained models, I use DeepHumanCortex1 here (snRNA+snATAC of 13 brain cell types).



