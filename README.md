# **Naïra Sarkis/Master Thesis**
In the following, further information on implemented tools and commands used for the bioinformatic analyses in my Master's Thesis are provided. All bioinformatic analyses have been performed on the HPC RRZK CHEOPS of the Regional Computing Centre (RRZK) of the University of Cologne using SLURM scripts, if not noted otherwise.

Versions of implemented programs:

![programs](https://user-images.githubusercontent.com/104494962/179422286-3e890e79-f4ad-4c1f-945c-b477e99df2fb.png)


## 1. *P. sambesii* whole genome annotation
### 1.1 Obtaining RNA-seq data from NCBI and preparing read files
> 
> The NCBI-SRA-toolkit is used to download paired-end Illumina HiSeq RNA-Seq data of *Plectus sambesii* from NCBI's SRA (SRR8243961). Note: This was done on a local computer and SRR8243961.tar.gz was scopied to CHEOPS.
> 
> ```prefetch  SRR8243961```
> 
> ```fasterq-dump -S SRR8243961``` 
> 
> For better program compatibility, headers of read files are modified using sed command. Dots in headers are exchanged with underspace and everything starting from the first blank space in a line is deleted.
> 
> ```tar -zcvf SRR8243961.tar.gz```
> 
> ```sed 's/\s.*$//' SRR8243961_1.fastq > SRR8243961_1.sed.fastq```
> 
> ```sed 's/\./_/g' SRR8243961_1.sed.fastq > SRR8243961_1.sednew.fastq```
> 
### 1.2 Masking the genome using repeatmodeler and repeatmasker
>
> ```BuildDatabase -name ES601_gene_DB -engine ncbi psambesii_genome.fasta RepeatModeler -engine ncbi -pa 16 -database ES601_gene_DB```
>
> ```RepeatClassifier -consensi ES601_gene_DB-families.fa```
>
> ```RepeatMasker -pa 16 -e ncbi -lib ES601_gene_DB-families.fa.classified psambesii_genome.fasta```
> 
### 1.3 Indexing genome and aligning reads using gmap/gsnap
>
> ```gmap_build -D /scratch/nsarkisk/Psam_annotation -d genome_index psambesii_genome.fasta.masked```
> 
> ```gsnap -D /scratch/nsarkisk/Psam_annotation -d genome_index -A sam -o /scratch/nsarkisk/Psam_annotation/psambesii-gsnap.sam SRR8243961_1.sednew.fastq SRR8243961_2.sednew.fastq```
>
### 1.4 Gene predictions using the braker2 pipeline
>
> ```braker.pl --species=PlectusSambesii --softmasking --AUGUSTUS_CONFIG_PATH=/scratch/nsarkisk/Psam_annotation/augustus-config/ --genome=psambesii_genome.fasta.masked --bam=psambesii-gsnap.bam.sorted```
> 
> A gff3 file containing the annotation is obtained and a file containing all coding sequences (CDS).

## 2. *P. sambesii* functional annotation 
### 2.1 Orthology Inference and hox gene analysis

### 2.1.1 OrthoFinder Analysis
>
> This section describes how the *P. sambesii* proteome was prepared for orthology analysis. 
>
> The annotation gff3 file from the braker output is converted into the right gff3 format using agat.
>
> ```agat_convert_sp_gxf2gxf.pl --gff psam_PB3_r3.braker.gff3 -o annotation.gff3```
>
> The longest peptide sequence is extracted.
> 
> ```agat_sp_extract_sequences.pl -g longest.gff3 -f psam-genome.fasta -o longest.fa -p```
> 
> Headers are changed to "PLESAM|ID" using sed command.
> 
> Directory "Fasta_files" is created and contains *Plectus sambesii* proteome from braker2 output and proteosome fasta files from species described in Material & Methods.
> 
> ```orthofinder -f Fasta_files/```
> 
### 2.1.2 Search orthogroups for hox genes
>
> The braker CDS ouput file had been translated using the python code python -m jcvi.formats.fasta beforehand to obtain a proteome file. A databank is generated from it using blast+.
> 
> ```makeblastdb -in psam_PB3_r3.braker3.fasta -dbtype prot -title Plectus-proteome```
> 
> Hox gene sequences (provided by Dr. Philipp Schiffer) were blasted against the new proteome.
> 
> ```blastp -query hox-proteins-plectus.fasta -db psam_PB3_r3.braker3.fasta -evalue 1e-30 -max_target_seqs 5 -outfmt 6 -out blastp_hox_vs_proteome.csv```
> 
> The output csv file contains hox genes as query and corresponding target sequence IDs in new *P. sambesii* proteome. Target sequence IDs are then used to find hox proteins in contigs of new annotation gff3 file.
>
> ```grep '>target-sequenceID<' psam_PB3_r3.braker.gff3 > hox_'target-sequenceID'```
>
> Annotation IDs of hox genes are used to search in output file "Orthogroups.txt" and saved as new files.
> Expl.: ```grep "g6951.t1" Orthogroups.txt > hox-g6951.t1```
> 
> OrthoFinder creates a directory "Orthogroup_Sequences" by default, containing amino acid fasta sequences of all genes within each orthogroup. Fasta files matching respective Hox IDs are stored in new directories.
> 
> Orthogroup fasta files to each hox protein are aligned using mafft.
> 
> ```mafft --localpair --maxiterate 1000 OG0000191.fa > OG0000191.fa.aln```
> 
> Spurious sequences are removed using trimal.
> 
> ```trimal -in OG0007202.fa.aln -out OG0007202.fa.aln.less.clw -resoverlap 0.75 -seqoverlap 80 -clustal```
> 
> Regions that do not align well are removed automatically.
>
> ```trimal -in OG0000191b.fa.aln -out OG0000191b.fa.aln.clean.clw -automated1 -clustal```
> 
> Maximum likelihood phylogenic trees are generated using iqtree.
> 
> ```iqtree2 -s OG0016393.fa.aln.clean.clw -m TEST -bb 1000 -nt AUTO```

### 2.2 Interproscan analysis

> The *P. sambesii* proteome file (from braker2) is split into 37 files using faSplit in order to run InterProscan with smaller files.
> 
> ```faSplit sequence psam_PB3_r3.braker3.aa 37 Psam```
> 
> InterproScan analysis is run for all chunked files.
> 
> ```singularity exec -B /home/nsarkisk/Psam_Interproscan_input/ -B /scratch/nsarkisk/Psam_interproscan_output/ -B /scratch/nsarkisk/interproscan-tmp/ /opt/rrzk/software/singularity_images/interproscan_latest.sif interproscan.sh -dp -T /scratch/nsarkisk/interproscan-tmp/ -goterms -i /home/nsarkisk/Psam_Interproscan_input/Psam01.fa -d /scratch/nsarkisk/Psam_interproscan_output/```
> 
> Chunked Interproscan results were combined to one tsv file.
> 
> ```cat Psam*.fa.tsv > interproscan-final.fa.tsv```

## 3. CEL-Seq2 pipeline
>
> This pipeline used process the sequencing results from CEL-Seq2 was developed by Dr. Tarja Hoffmeyer of the worm~lab, Cologne.
>
### 3.1. Removing adapters and quality control
>
> Create seperate fasta files with known adapters that can be looked up in the fastp known adapters file (https://github.com/OpenGene/fastp/blob/master/src/knownadapters.h).
> 
> Create seperate fasta files for each Library in this format:
>
> ![Psam_adapters](https://user-images.githubusercontent.com/104494962/177975386-41e916b1-57fe-4c07-ad79-82997981312c.png) 
> 
> RP1 is the universal adapter for small RNA Illumina libraries and RPI2, RP9, RPI10, RPI11 are the indexed adapters, respectively.
>
> Adapters are removed according to defined adapters in fasta files, followed by quality control.
> 
> ```fastp -i Psam-1_1.fq.gz -I Psam-1_2.fq.gz -o Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_forward_paired.fq.gz -O Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_reverse_paired.fq.gz --unpaired1 Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_forward_unpaired.fq.gz --unpaired2 Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_reverse_unpaired.fq.gz -g --poly_g_min_len 3 -x -l 12 --adapter_fasta adapters_Psam-1.fa --cut_front --cut_front_mean_quality 3 --cut_tail --cut_tail_mean_quality 3 --cut_right --cut_right_mean_quality 15 -p -j Psam-1_fastp_adapter_fasta_polyg3_polyx_min12.json -h Psam-1_fastp_adapter_fasta_polyg3_polyx_min12.html```
>
> In NextSeq, a polyG tail often occurs at the end of reads as unreadable bases are read as G. -g flag these if least three Gs in a row are present. Other poly tails are removed with -x option. Reads with less than 12 bases are removed with option -l 12. Bad quality front and tail areas are removed. The option --cut_right is a sliding window function that removes very low quality regions in the middle areas. Overrpresented read analysis in performed with -p flag. -json and -html outputs are created. A minimum length of 12 is given for trimming to remove reads only containing UMI and barcode.
>
> Reads shorter than 36 bases are removed from reverse paired reads.
> 
> ```fastp -i Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_reverse_paired.fq.gz -o Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_reverse_paired.fq.gz -l 36 -p -j Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_reverse_paired.json -h Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_reverse_paired.html```
>
> Removing low complexity reads (default complexity level).  
>
> ```fastp -i Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_reverse_paired.fq.gz -o Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_reverse_paired.fq.gz -l 36 -y -p -j Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_reverse_paired.json -h Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_reverse_paired.html```
>
> Removing reads with 15xA bases using bbduk.
>
> ```bbduk.sh k=15 in=Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_reverse_paired.fq.gz out=Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_bbduk_minuspolyA_k15_hdist0_reverse_paired.fq.gz outm=Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_bbduk_onlypolyA_k15_hdist0_reverse_paired.fq.gz literal=AAAAAAAAAAAAAAA hammingdistance=0```
>
> Remaining reverse reads are paired with forward reads and forward reads without pair are removed. A file with the fraction of the headers of the reverse paired file that can also be found in the corresponding headers of the forward read file is prepared for all libraries.
>
> ```for f in Psam-1 Psam-3 Psam-4 Psam-5; do sed -n '1~4p' $f.fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_bbduk_minuspolyA_k15_hdist0_reverse-paired.fq | sed 's/\s.*$//' | sed 's/^@//g' > $f.final_bbduk_reverse_paired_headers; done```
>
> Sequences from reverse paired file of the forward paired file are saved into a new file using pullseq. 
>
> ```for f in Psam-1 Psam-3 Psam-4 Psam-5; do pullseq -i ../$f.fastp_adapter_fasta_polyg3_polyx_min12_forward_paired.fq -n $f.final_bbduk_reverse_paired_headers > $f.final_bbduk_forward_paired.fq; done```
>
> Repeat for reverse reads.
>
> ```for f in Psam-1 Psam-3 Psam-4 Psam-5; do pullseq -i ../$f.fastp_adapter_fasta_polyg3_polyx_min12_reverse_paired.fq -n $f.final_bbduk_reverse_paired_headers > $f.final_bbduk_reverse_paired.fq; done```
>
### 3.2. Extract UMIs from all forward reads and write into header of both forward and reverse reads
>
> ```for f in Psam-1 Psam-3 Psam-4 Psam-5; do umi_tools extract -I $f.finaltrim_forward_paired.fq --bc-pattern=NNNNNN --read2-in=$f.finaltrim_reverse_paired.fq --stdout=$f.finaltrim_forward_paired_umiext.fq --read2-out=$f.finaltrim_reverse_paired_umiext.fq; done```
>
> Run demultiplex to sort the 12 bases forward fastq reads into separate files according to their Celseq2 barcode. Bases are sorted at positions 0-6 (UMIs have been extracted), while the programme allows for one mismatch by default. The barcodes are all different in at least two positions and just one error will still allow to sort the reas to the respective barcode.
> 
> The 24 barcodes that were used are written into an Excel file and saved as tab delimited .txt file and renamed with suffix .tsv [Barcodes](https://user-images.githubusercontent.com/104494962/177957129-58f7e0f9-d799-4981-88a7-6e7a1bfc2f9a.png)
> 
> The forward paired reads after fastp trimming are used and sorted into files according to the embryo sample they came from by using the barcodes.
> 
> ```for f in Psam-1 Psam-3 Psam-4 Psam-5; do demultiplex demux -r -s 1 -e 6 barcodes24.tsv $f.finaltrim_forward_paired_umiext.fq; done```
> 
> Only the headers from the forward reads are needed and are sorted into files according to the embryo sample they came from using sed.
> 
> ```for f in Psam-1.finaltrim_forward_paired_Celseq*.fq; do sed -n '1~4p' $f > $f.list; done``` Repeat for each library.
>
> Create a list of all headers without "@" and without everything after the first space character.
> 
>```for f in *.list; do sed 's/\s.*$//' $f > $f.fixedlist; done```
>
>```for f in *.list.fixedlist; do sed 's/^@//g' $f > $f.2; done```
>
> The sequences for each sample are extracted into seperate files.
> 
> ```for f in Psam-1*.list.fixedlist.2; do pullseq -i Psam-1_fastp_reverse_paired_umiext.fq -n $f > $f.sequences_pullseq.fq; done``` Repeat for each library.
>
### 3.3. Mapping reads onto genome
>
> Indexing the genome using kallisto. k-mer size is set 21 for small reads. 
> 
> ```kallisto index -i P_sambesii21.index -k 21 psambesii_genome.fasta.masked```
> 
> Mapping of each of the demultiplexed samples against the index. Only transcript read is mapped. -l flag describes size of library (200bp, as measured by femto pulse analysis), -s flag describes variation, set to 10 % of the size here. Pseudomapping is performed to retain the information of how many reads map to which gene.
> 
> ```for f in ./Psam-1*sequences_pullseq.fq; do kallisto quant -i ./P_sambesii21.index --pseudobam -o $f.kallisto21 --single -l 200 -s 20 $f > $f.kallisto21.sam; done; for f in ./Psam-3*sequences_pullseq.fq; do kallisto quant -i ./P_sambesii21.index --pseudobam -o $f.kallisto21 --single -l 200 -s 20 $f > $f.kallisto21.sam; done; for f in ./Psam-4*sequences_pullseq.fq; do kallisto quant -i ./P_sambesii21.index --pseudobam -o $f.kallisto21 --single -l 200 -s 20 $f > $f.kallisto21.sam; done; for f in ./Psam-5*sequences_pullseq.fq; do kallisto quant -i ./P_sambesii21.index --pseudobam -o $f.kallisto21 --single -l 200 -s 20 $f > $f.kallisto21.sam; done```
> 
### 3.4. Sorting and indexing BAM files
> 
> Since all created kallisto subdirectories are named "pseudoalignments.bam", they are renamed and retrieved from the subdirectories. All mapping information is stored in the pseudobam files.
> 
> ```for f in *.bam; do samtools view -b $f > $f.view; done```
> ```for f in *.view; do samtools sort $f -o $f.sorted.bam; done```
> ```for f in *.sorted.bam; do samtools index $f; done```
> 
> The created .bai index files need to be in the same directory moving forward.
>
> The mapping statistics were summarized by using MultiQC and retrieved as .html file.
>
> ```multiqc .```
>
### 3.5. UMI-tools deduplication
> 
> Create a list of files
> 
> ```find ./*kallisto21.bam.view.sorted.bam > kallisto21_sorted_bam_list```
> 
> Deduplication using umi_tools
> 
> ```cat kallisto21_sorted_bam_list | parallel -j 12 umi_tools dedup -I {} --output-stats={}.deduplicated -S {}.deduplicated.bam```
> 
### 3.6. Create expression matrix from UMI counts
> 
> Create a .gtf file from the .gff3 annotation file
> 
> ```agat_convert_sp_gff2gtf.pl --gff psam_PB3_r3.braker.gff3 -o reference_genome_annotation.gtf```
> 
> Assign reads to genes using featureCounts from Subread package
> 
> ```for f in *deduplicated.bam; do featureCounts -a reference_genome_annotation.gtf -o $f.genes_assigned -R BAM $f -T 4; done```
> 
> Use samtools to sort and index
> 
> ```for f in *.bam.featureCounts.bam; do samtools sort $f -o $f.assigned_sorted.bam; done```
> 
> ```for f in *.assigned_sorted.bam; do samtools index $f; done```
> 
> Count genes using umi_tools and unzip all
> 
> ```for f in *.assigned_sorted.bam; do umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I $f -S $f.counts.tsv.gz; done```
> 
> ```gunzip *.counts.tsv.gz```
> 
> Counts were made for cells, but each sample was for one embryo, so the full counts per embryo are needed. Start with only extracting the first line with the repeating gene names
> 
> ```for f in *.counts.tsv; do awk '{print $1}' $f > $f.genes.tsv; done```
> 
> Sort the gene names, only keep the unique ones and count the occurrences
>
> ```for f in *.counts.tsv.genes.tsv; do sort $f | uniq -c > $f.sorted_counted.tsv; done```
> 
> Files have counts in the first column and genes names in the second column. Swap the two columns
> 
> ```for f in *.sorted_counted.tsv; do awk '{t=$1; $1 =$2; $2=t;print;}' $f > $f.swapped.tsv; done```
> 
> Create new directory to store all count files and replace all spaced in the files with tabs
> 
> ```mkdir counts```
> 
> ```for f in *.swapped.tsv; do sed 's/\s/\t/g' $f > ./counts/$f.tabs.tsv; done```
> 
> Rename all files in the directory to this form: Lib-1.Celseq1.tsv
>
> ```rename 's/_fastp_forward_paired_umiext_/./' *.tsv```
> 
> ```rename 's/.fq.list.fixedlist.2.sequences_pullseq.fq.kallisto21.bam.view.sorted.bam.deduplicated.bam.featureCounts.bam.assigned_sorted.bam.counts.tsv.genes.tsv.sorted_counted.tsv.swapped.tsv.tabs.tsv/.tsv/' *.tsv```
> 
> Delete all header lines
>
> ```for f in *.tsv; do sed -i 1d $f; done```
>
> Show and remove completely empty files and move only final counts to a new directory
>
> ```find . -type f -empty -print -delete```
>
> ```mv *.tsv ./new_directory/```
>
> Create a conda environment with R version 4.1 and install R packages
> 
> ```conda create -n r_4_1 -c conda-forge r-base=4.1* r-essentials```
> 
> ```conda activate r_4_1```
> 
> ```R```
> 
> ```if (!require("BiocManager", quietly =TRUE)) install.packages("BiocManager")```
> 
> ```BiocManager::install("SingleCellExperiment")```
> 
> ```BiocManager::install("scuttle")```
> 
> ```BiocManager::install("scater")```
> 
> ```BiocManager::install("edgeR")```
> 
> ```install.packages("Matrix")```
> 
> ```install.packages("dplyr")```
>
> Load required packages
>
> ```library(dplyr)```
> 
> ```library(Matrix)```
> 
> ```library(SingleCellExperiment)```
> 
> ```library(scuttle)```
> 
> ```library(scater)```
> 
> ```library(edgeR)```
>
> Create a character list of prepared UMI count files from within the directory using R
>
> ```files <- list.files()``` Creates vector from list of file names in the working directory
> 
> ```UMI_counts <- readDGE(files)``` Reads and stores the data from all files stored in file vector in UMI_counts
>
> ```ountData <- as.data.frame(UMI_counts$counts)``` Creates count matrix
> 
> ```write.csv(countData, "./UMI_counts.csv", row.names = TRUE)``` Write count matrix to .csv file
>
> ```data = read.delim("./UMI_counts.csv", sep=",", header=TRUE)``` Read in the data from the file as a data frame
> 
> ```data[1:5,1:5]``` Check first 5 rows and columns of the dataframe
> 
> ```dim(data)``` Check dimensions of the data frame
> 
> ```rownames(data) <- data[,1]``` Set rownames to the first column
> 
> ```countsmatrix <- data[,-1]``` Remove the first column, as the names are now stored in rownames(data). The dataframe is now numeric and can be transformed into a matrix
> 
> ```options(stringsAsFactors = FALSE)``` Global default setting changed; every dataframe created hereafter will not auto-convert to factors unless explicitly told to do so
> 
> ```umi <- SingleCellExperiment(assays=list(counts=as.matrix(countsmatrix)))``` Creates sce object (S4 class for storing data from single-cell experiments. This object can be used by the package scater for conversion, QC, and normalisation
> 
> ```umi <- umi[, colSums(counts(umi)) > 0]``` Spots with non-positive counts are removed
> 
> ```tpm(umi) <- calculateTPM(umi, length=NULL, assay.type="counts", exprs_values=NULL)``` Conversion of counts to TPM (scuttle function), length set to NULL for UMI counts
> 
> ```tpm(umi)[1:5,1:5]``` Check first 5 rows and columns
> 
> ```normcounts(umi) = log10(tpm(umi)+1)``` log10 normalisation (+1 added to not have null values; cannot be log transformed)
> 
> ```write.csv(normcounts(umi), "./UMI_tpm_log10.csv", row.names = TRUE)``` Save as UMI_tpm_log10.csv file
> 
> ```keep_feature <- rowSums(normcounts(umi) > 0) > 0```
> 
> ```umi2 <- normcounts(umi)[keep_feature, ]``` All entries with just zeros are removed
> 
> ```UMIcountsperembryo <- colSums(umi2)``` Sums together all counts in each column
> 
> ```genecountsperembryo <- colSums(umi2 !=0)``` Sums together all fields in each column that are not 0
> 
> ```nGenesrationUMI <- genecountsperembryo/UMIcountsperembryo``` Calculates ratio
> 
> ```write.csv(UMIcountsperembryo, "./UMIcountsperembryo.csv", row.names = TRUE)``` Writes .csv file
> 
### 3.7. QC of CEL-Seq2 single embryo experiment
>  
> UMI counts per embryo, gene counts per embryo and complexity were determined in the last step. Here the ratio between UMIs corresponding to mitochondrial genes vs. other genes is calculated. This is done more easily if the mitochondrial contig is part of the assembly. In this thesis an genome assembly was used, where the mitochondrial contig had been removed in a decontamination step. The mitochondrial contig is being retrieved from an original assembly via BLAST search against a database created for this assembly.
> 
> Create a local database for the raw genome assembly
> ```module add blast+```
> ```makeblastdb -in raw_genome_assembly.fasta -dbtype nucl -title raw_genome_assembly```
> 
> BLAST search of the C. elegans cox-1 protein to identify the mitochondrial contig.
> 
> ```tblastn -query CELE_cox1aa.fa -db raw_genome_assembly.fasta -evalue 1e-2 -max_target_seqs 5 -outfmt 6 -out```
> 
> Extract the contig and remove lines inbetween the sequences.
> 
> ```awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' raw_genome_assembly.fasta > raw_genome_assembly_minus_new_lines.fasta```
> 
> Check for the line containing the string ">contig_3238" and copy this and the next line containing the sequence into a new file.
> 
> ```grep -A 1 ">contig_3238" raw_genome_assembly_minus_new_lines.fasta > contig3238.fasta```
> 
> Confirm the identity of the mitochondrial contig by perfomring a BLAST search against the NCBI nt database.
> 
> MITOS (http://mitos.bioinf.uni-leipzig.de/) is used to annotate the mitochondrial contig. An .gff annotation file is created and shows all regions present in the mitochondrial contig.
> 
> In order to retrieve the mitochondrial UMI counts, mapping and featureCounts is repeated as described before for the mitochondrial contig. For featureCounts of the mitochondrial counts the -t and -g flag has to be added.
> 
> ```for f in *deduplicated.sorted.bam; do featureCounts -a Consensusmito.gtf -o $f.genes_assigned -t gene -g gene_id -R BAM $f -T 4; done```
> 
> Create a table in Microsoft Excel in this format:
> Embryo   nUMI  nGenes  nUMI/nGenes  mito UMI fraction
> 
### 3.8. Sorting of gene expression and heatmap creation
> 
> Same R environment is used as above.
> 
> ```conda activate r_4_1```
> 
> ```R```
> 
> ```library(dplyr)```
> 
> ```library(Matrix)```
> 
> ```library(SingleCellExperiment)```
> 
> ```library(scuttle)```
>
> ```library(scater)```
>
> ```library(edgeR)```
> 
> ```UMI_tpm_log10 = read.delim("./UMI_tpm_log10.csv", sep=",", header=TRUE)```
> 
> ```rownames(UMI_tpm_log10) <- UMI_tpm_log10[,1]``` Determine rownames
> 
> ```UMI_tpm_log10_counts <- UMI_tpm_log10[,-1]``` Remove column containing rownames as they are stored anyway
> 
> ```col_order = c("Psam.3_Celseq1", "Psam.1_Celseq17", "Psam.1_Celseq1", "Psam.3_Celseq80", "Psam.3_Celseq31", "Psam.5_Celseq26", "Psam.1_Celseq80", "Psam.4_Celseq24", "Psam.3_Celseq77", "Psam.4_Celseq74", "Psam.1_Celseq50", "Psam.4_Celseq44", "Psam.3_Celseq63", "Psam.3_Celseq74", "Psam.1_Celseq68", "Psam.5_Celseq17", "Psam.5_Celseq55", "Psam.5_Celseq23", "Psam.1_Celseq26", "Psam.4_Celseq96", "Psam.4_Celseq10", "Psam.4_Celseq9", "Psam.1_Celseq5", "Psam.1_Celseq6", "Psam.1_Celseq44", "Psam.5_Celseq1", "Psam.3_Celseq6", "Psam.1_Celseq31", "Psam.5_Celseq80", "Psam.5_Celseq24", "Psam.4_Celseq1", "Psam.4_Celseq36", "Psam.5_Celseq63", "Psam.3_Celseq17", "Psam.4_Celseq89", "Psam.1_Celseq89", "Psam.4_Celseq4", "Psam.5_Celseq4", "Psam.4_Celseq6", "Psam.3_Celseq4", "Psam.3_Celseq9", "Psam.5_Celseq68", "Psam.3_Celseq68", "Psam.4_Celseq50", "Psam.1_Celseq10", "Psam.1_Celseq23", "Psam.1_Celseq55", "Psam.4_Celseq46", "Psam.5_Celseq36", "Psam.5_Celseq10", "Psam.3_Celseq89", "Psam.3_Celseq50", "Psam.5_Celseq74", "Psam.3_Celseq10", "Psam.3_Celseq25", "Psam.5_Celseq50", "Psam.5_Celseq89", "Psam.3_Celseq26", "Psam.4_Celseq68", "Psam.3_Celseq23", "Psam.3_Celseq5", "Psam.1_Celseq36", "Psam.1_Celseq74", "Psam.4_Celseq26", "Psam.3_Celseq55", "Psam.4_Celseq5", "Psam.5_Celseq46", "Psam.1_Celseq24", "Psam.5_Celseq25", "Psam.5_Celseq77", "Psam.3_Celseq24", "Psam.4_Celseq80", "Psam.4_Celseq77", "Psam.4_Celseq17", "Psam.1_Celseq25", "Psam.5_Celseq44", "Psam.3_Celseq44", "Psam.1_Celseq4", "Psam.1_Celseq77", "Psam.1_Celseq63", "Psam.3_Celseq36", "Psam.1_Celseq46", "Psam.1_Celseq9", "Psam.5_Celseq5", "Psam.3_Celseq46", "Psam.4_Celseq23", "Psam.5_Celseq96", "Psam.5_Celseq6", "Psam.5_Celseq31", "Psam.5_Celseq9", "Psam.4_Celseq25", "Psam.3_Celseq96", "Psam.4_Celseq63", "Psam.4_Celseq31", "Psam.4_Celseq55")```
> 
> ```UMI_tpm_log10_ordered <- UMI_tpm_log10_counts[,col_order]``` Reorders the columns of the entire table
> 
> ```write.csv(UMI_tpm_log10_ordered, "./UMI_tpm_log10_ordered.csv", row.names = TRUE)``` Save file with new order
> 
> ```install.packages("rgdal")```
> 
> ```install.packages("iemisc")```
> 
> Load the normalized expression matrix
> 
> ```normexp = read.delim("./UMI_tpm_log10_ordered1.csv", sep=",", header=TRUE)```
> 
> ```rownames(normexp) <- normexp[,1]``` Creates a file containing the rownames of the data frame
> 
> ```normexpmatrix <- normexp[,-1]``` Remove rownames from the dataframe to be able to convert it to a matrix
> 
> ```mat <- as.matrix(normexpmatrix)``` Convert to a matrix
> 
> Zavit code (received from Gustavo Starvaggi Franco) is used with this matrix in R
> 
> ```exp_tpm_scale = t(scale(t(mat)))```
> ```exp_tpm_pca = princomp(exp_tpm_scale, cor=F)```
> ```X = exp_tpm_pca$scores[, c("Comp.1")]```
> ```Y = exp_tpm_pca$scores[, c("Comp.2")]```
> ```library("iemisc")```
> ```t = atan2d(X, Y)```
> ```t_ordered = sort(t, decreasing = T)```
> ```exp_tpm_ord = mat[names(t_ordered), ]```
> 
> ```write.csv(exp_tpm_ord, "./exp_tpm_ord.csv", row.names = TRUE)``` Creates a file with the new gene order
> 
> Create a heatmap from the same R session
> 
> ```library(ggplot2)```
> 
> ```library(reshape2)```
> 
> ```df <- read.csv("exp_tpm_ord.csv")```
> 
> ```dfmelt <- melt(df)``` Restructures the data, so that ggplot can use it
> 
> ```heatmapdf <- ggplot(dfmelt, aes(x = variable, y = X, fill = value)) + geom_tile() + scale_fill_gradient(high = "red", low = "blue")``` Creates heatmap
> 
> ```pdf("heatmapdf.pdf")``` Creates a PDF file
> 
> ```print(heatmapdf)``` Prints the heatmap to the PDF file
> 
> ```dev.off()``` Closes the PDF file
