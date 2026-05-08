# TF Analysis
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
Make use of CRESted and TF-MINDI for the prediciton of motif syntax disruption by variant in a cell-type-specific manner.

### Setup

Pip install TF-MINDI and CRESted[motif] (with Tensorflow and cuda), modisco (motif discovery)
Dependendies: pybedtools pysam

There's two human brain pre-trained models, I use DeepHumanCortex1 here (snRNA+snATAC of 13 brain cell types). Downloand and check with brain_model_crested.py

### CRESted prediction for open chromatin from sequence alone
Run crested_contr_score.py for calculating the importance of each nulceotide for the prediction of chromatin accessibility in a cell type.
1. Extract ±250bp reference and alternative sequences from FASTA (model was trained on 500bp so cannot get higher).
2. One-hot encode both
3. Compute per-nucleotide contribution scores (once per sequence)
4. Save per-SNP directory with one-hot arrays and contribution scores
   
If calculated, the delta between REF-ALT nucleotide contribution score would mean:
Positive delta -> change in the allele boosted the model's prediction for open chromatin in this region 
Negative delta -> change in the allele weaked the model's prediction for open chromatin in this region

### Plotting contribution scores 
The contribution scores (per cell type) can be visualized with alt_vs_ref_plot.py. This will produce two figures, one reference and one alternate sequence plot with SNP in the centre. 

### Run TF-MINDI motif similarity analysis
Run TF_analysis_across_all.py. 
TF-MINDI extracts seqlets (important sequence fragments) based on the contribution score. Here, however, I am extracting the 20bp window around my SNP of interest, whether it is a seqlet or not. Motif similarity will scan the sequence syntax for familiar motifs levaraging the cell-type-specific contribution scores of neighboring nucleotides. This gives the analysis context-dependent and cell-type-specific advantage. From a strip of sequence, different motifs can match to different cell types, based on their dissimilar contribution scores. Delta is calculated from REF-ALR motif similarity and further statistical corrections are added. The script categorizes variants as
1. Gain -> motif syntax strengthened upon alternate allele
2. Loss -> motif syntax weakedn upon alternate allele

Summarize results by running tf_family_summary.py
Plot the top motif gain/loss in a table with plot_per_snp_tf_family_table.py
If you have a managable list of SNPs of interest, you can make a nice heatmap to see gain/loss motifs in every cell type by running plot_snp_celltype_tf_heatmap.py











