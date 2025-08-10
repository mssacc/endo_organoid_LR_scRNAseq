###FLAMES analysis of all samples

#Load Required Libraries
library("FLAMES")

#Config File
config_file <- "/data/gpfs/projects/punim1901/flames_v2/config.json"

#Fastq Files
fastqs = c("/data/gpfs/projects/punim1901/Endo_organoid_multisample/fastq/E170/E170.fastq.gz",
           "/data/gpfs/projects/punim1901/Endo_organoid_multisample/fastq/E191/E191.fastq.gz", 
           "/data/gpfs/projects/punim1901/Endo_organoid_multisample/fastq/E226/E226.fastq.gz",
           "/data/gpfs/projects/punim1901/Endo_organoid_multisample/fastq/E231/E231.fastq.gz", 
           "/data/gpfs/projects/punim1901/Endo_organoid_multisample/fastq/E333/E333.fastq.gz", 
           "/data/gpfs/projects/punim1901/Endo_organoid_multisample/fastq/E435/E435.fastq.gz")

#Output Directory
output = "/data/gpfs/projects/punim1901/flames_v2/flames_output"

#Reference Files
GTF = "/data/gpfs/projects/punim1901/genomes/gencode.v31.annotation.gtf"
genome = "/data/gpfs/projects/punim1901/genomes/hg38.analysisSet.fa"
	
sce <- sc_long_multisample_pipeline(fastqs=fastqs, outdir=output, annot=GTF, genome_fa=genome, 
config_file=config_file, expect_cell_numbers=c(2000,2000,2000,2000,2000,2000))
