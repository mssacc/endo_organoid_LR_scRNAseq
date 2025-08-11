# Endometrial epithelial organoid long-read single-cell RNA-sequencing
This repository contains scripts required to produce results and figures presented in paper "Long-read single-cell RNA sequencing of endometrial organoids uncovers cell-type specific isoform dysregulation in unexplained infertility"

Long-read scRNA-seq fast5, pod5, and fastq files are publicly available for download from ENA under: (add id)

General data analysis outline
1) Run FLAMES multisample analysis
2) Run FLAMES empty droplet analysis
3) Generate reference files
4) Convert gene ID to gene symbol
5) Run gene QC filtering
6) Combine gene and isoform data
7) Merge and integrate samples

Additional analyses
1) Cell subtype composition and variability
2) Cell-cell communication
3) Differential gene and isoform expression (global + cell subtype specific)
4) Isoform classification
5) Differential transcript isoform usage
6) Transcription factor assessment
7) Gene set enrichment analysis
