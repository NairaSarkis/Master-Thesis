# **Naïra Sarkis/Master Thesis**
In the following, further information on implemented tools and commands used for the bioinformatic analyses in my Master's Thesis are provided. All bioinformatic analyses have been performed on the HPC RRZK CHEOPS of the Regional Computing Centre (RRZK) of the University of Cologne, if not noted otherwise. 

Versions of implemented programs:
###Table
## 1. *P. sambesii* whole genome annotation
### 1.1 Obtaining RNA-seq data from NCBI and preparing read files
> The NCBI-SRA-toolkit was used to download paired-end RNA-seq data for *Plectus sambesii* that has been generated before (###Ref)
> 
> ```fasterq-dump -S SRR8243961```
> 
> For better program compatibility, headers of read files were modified using sed command. Dots in headers were exchanged with underspace and everything starting from the first blank space in a line was deleted.
> 
### 1.2 Masking the genome using repeatmodeler and repeatmasker
> ```BuildDatabase -name ES601_gene_DB -engine ncbi psambesii_genome.fasta RepeatModeler -engine ncbi -pa 16 -database ES601_gene_DB```
>
> ```RepeatClassifier -consensi ES601_gene_DB-families.fa```
>
> ```RepeatMasker -pa 16 -e ncbi -lib ES601_gene_DB-families.fa.classified psambesii_genome.fasta```
### 1.3 Indexing genome and aligning reads using gmap/gsnap
> ```gmap_build -D /scratch/nsarkisk/Psam_annotation -d genome_index psambesii_genome.fasta.masked```
> 
> ```gsnap -D /scratch/nsarkisk/Psam_annotation -d genome_index -A sam -o /scratch/nsarkisk/Psam_annotation/psambesii-gsnap.sam SRR8243961_1.sednew.fastq SRR8243961_2.sednew.fastq```
### 1.4 Gene predictions using the braker2 pipeline
> ```braker.pl --species=PlectusSambesii --softmasking --AUGUSTUS_CONFIG_PATH=/scratch/nsarkisk/Psam_annotation/augustus-config/ --genome=psambesii_genome.fasta.masked --bam=psambesii-gsnap.bam.sorted```
> 
> A gff3 file containing the annotation was obtained and a file containing all coding sequences (CDS).

## 2. *P. sambesii* functional annotation 
### 2.1 Orthology Inference and hox gene analysis

2.1.1 OrthoFinder Analysis
> 
> The annotation gff3 file from the braker output is converted into the right gff3 format using agat
> 
> ```agat_sp_extract_sequences.pl -g longest.gff3 -f psam-genome_folded.fasta -o longest.fa -p```
> 
> Headers were changed to "PLESAM|ID" using sed command.
> Directory "Fasta_files" is created and contains proteosome fasta files from ###SPECIES and Plectus sambesii proteome from braker2 output.
> 
> ```orthofinder -f Fasta_files/```
> 
2.1.2 Search orthogroups for hox genes
> The braker CDS ouput file is translated in order to get a proteome file (###method) and a databank is generated from it using blast+. 
> 
> ```makeblastdb -in psam_PB3_r3.braker3.fasta -dbtype prot -title Plectus-proteome```
> Previously provided hox gene sequences (###ref) were blasted against the new proteome.
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

> The P. sambesii proteome file (from braker2) was split into 37 files using faSplit in order to run InterProscan with smaller files.
> 
> ```faSplit sequence psam_PB3_r3.braker3.aa 37 Psam```
> 
> ```for f in Psam*.fa; do /my_interproscan/interproscan-5.55-88.0/interproscan.sh -i $f -f tsv --goterms --pathways; done```
> 
> Chunked Interproscan results were combined to one tsv file.
> 
> ```cat Psam*.fa.tsv > interproscan-final.fa.tsv```

## 3. CELSeq2 downstream analysis 

### Tools implemented in this analysis:

> ###Table of versions
### 1. Removing adapters and quality control
> 1.1. Create fasta files with known adapters
> 
> 1.2. Adapters are removed according to defined adapters in fasta files, followed by quality control
> 
> ```fastp -i Psam-1_1.fq.gz -I Psam-1_2.fq.gz -o Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_forward_paired.fq.gz -O Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_reverse_paired.fq.gz --unpaired1 Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_forward_unpaired.fq.gz --unpaired2 Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_reverse_unpaired.fq.gz -g --poly_g_min_len 3 -x -l 12 --adapter_fasta adapters_Psam-1.fa --cut_front --cut_front_mean_quality 3 --cut_tail --cut_tail_mean_quality 3 --cut_right --cut_right_mean_quality 15 -p -j Psam-1_fastp_adapter_fasta_polyg3_polyx_min12.json -h Psam-1_fastp_adapter_fasta_polyg3_polyx_min12.html```
>
> 1.3. Reads shorter than 36 bases are removed from reverse paired reads
> 
> ```fastp -i Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_reverse_paired.fq.gz -o Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_reverse_paired.fq.gz -l 36 -p -j Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_reverse_paired.json -h Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_reverse_paired.html```
>
> 1.4. Removing low complexity reads
>
> Complexity describes the percentage of base that is different from its next base. The complexity filter default value is 30, which means 30 % complexity is required, to keep the reads. (###SOURCE?)
>
> ```fastp -i Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_reverse_paired.fq.gz -o Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_reverse_paired.fq.gz -l 36 -y -p -j Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_reverse_paired.json -h Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_reverse_paired.html```
>
> 1.5. Removing reads with 15xA bases using bbduk
>
> ```bbduk.sh k=15 in=Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_reverse_paired.fq.gz out=Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_bbduk_minuspolyA_k15_hdist0_reverse_paired.fq.gz outm=Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_bbduk_onlypolyA_k15_hdist0_reverse_paired.fq.gz literal=AAAAAAAAAAAAAAA hammingdistance=0```
>
> 1.6. Remaining reverse reads are paired with forward reads again
>
> Forward reads that are missing are removed. 
>
> 1.7. Create a file, that contains a fraction of headers of the reverse paired file, that can also be found in the corresponding headers of forward read files
>
> ```zcat Psam-1_fastp_adapter_fasta_polyg3_polyx_min36_minlowcomplexity_bbduk_minuspolyA_k15_hdist0_reverse_paired.fq.gz | sed -n '1~4p' | sed 's/\s.*$//' | sed 's/^@//g' > Psam-1_final_bbduk_reverse_paired_headers```
>
> 1.8. Pull sequences from reverse paired file of the forward paired file into a new file using pullseq
>
> ```pullseq -i Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_forward_paired.fq -n Psam-1_final_bbduk_reverse_paired_headers > Psam-1_final_bbduk_forward_paired.fq```
>
> Repeat for reverse reads
>
> ```pullseq -i Psam-1_fastp_adapter_fasta_polyg3_polyx_min12_reverse_paired.fq -n Psam-1_final_bbduk_reverse_paired_headers > Psam-1_final_bbduk_reverse_paired.fq```
>
### 2. Extract UMIs from all forward reads and write into header of both forward and reverse reads
> ```umi_tools extract -I Psam-1_fastp_forward_paired.fq --bc-pattern=NNNNNN --read2-in=Psam-1_fastp_reverse_paired.fq --stdout=Psam-1_fastp_forward_paired_umiext.fq --read2-out=Psam-1_fastp_reverse_paired_umiext.fq```
>
> Run demultiplex to sort the 12 bases forward fastq reads into separate files according to their Celseq2 barcode.
> We are sorting the bases at positions 0-6 (UMIs have been extracted). The programme, by default, allows for one mismatch.
> 
> The 24 barcodes that were used are written into an Excel file and saved as tab delimited .txt file and renamed with suffix .tsv
> [Barcodes](https://user-images.githubusercontent.com/104494962/177957129-58f7e0f9-d799-4981-88a7-6e7a1bfc2f9a.png)

> The forward paired reads are used and sorted into files according to the embryo sample they came from by using the barcodes.
> 
> ```demultiplex demux -r -s 1 -e 6 barcodes24.tsv Psam-1_fastp_forward_paired_umiext.fq```
> 
> Only the headers from the forward reads are needed and extracted using sed.
> 
> ```for f in Psam-1_fastp_forward_paired_Celseq*.fq; do sed -n '1~4p' $f > $f.list; done```
>
> Create a list of all headers without "@" and without everything after the first space character.
> 
>```for f in *.list; do sed 's/\s.*$//' $f > $f.fixedlist; done```
>
>```for f in *.list.fixedlist; do sed 's/^@//g' $f > $f.2; done```
>
> The sequences for each sample are extracted into seperate files
> 
> ```for f in Psam-1*.list.fixedlist.2; do pullseq -i Psam-1_fastp_reverse_paired_umiext.fq -n $f > $f.sequences_pullseq.fq; done```
>
### 3. Mapping reads onto genome
>
> Indexing the genome using kallisto. k-mer size is set 21 for small reads. 
> 
> ```kallisto index -i P_sambesii21.index -k 21 psambesii_genome.fasta.masked```
> 
> Mapping of each of the demultiplexed samples against the index. Only transcript read is mapped. -l flag describes size of library (200bp, as measured by femto pulse analysis), -s flag describes variation, set to 10 % of the size here. 
> 
> ```for f in ./Psam-1*sequences_pullseq.fq; do kallisto quant -i ./P_sambesii21.index --pseudobam -o $f.kallisto21 --single -l 200 -s 20 $f > $f.kallisto21.sam; done; for f in ./Psam-3*sequences_pullseq.fq; do kallisto quant -i ./P_sambesii21.index --pseudobam -o $f.kallisto21 --single -l 200 -s 20 $f > $f.kallisto21.sam; done; for f in ./Psam-4*sequences_pullseq.fq; do kallisto quant -i ./P_sambesii21.index --pseudobam -o $f.kallisto21 --single -l 200 -s 20 $f > $f.kallisto21.sam; done; for f in ./Psam-5*sequences_pullseq.fq; do kallisto quant -i ./P_sambesii21.index --pseudobam -o $f.kallisto21 --single -l 200 -s 20 $f > $f.kallisto21.sam; done```
> 
> ### 4. Sorting and indexing BAM files
> 
> Since all created kallisto subdirectories are named "pseudoalignments.bam", they were renamed and retrieved from the subdirectories. All mapping information is stored in the pseudobam files.
> 
> ```for f in *.bam; do samtools view -b $f > $f.view; done```
> ```for f in *.view; do samtools sort $f -o $f.sorted.bam; done```
> ```for f in *.sorted.bam; do samtools index $f; done```
> 
> The created .bai index files need to be in the same directory moving forward.
>
> ### 5. UMI-tools deduplication
> 
> Create a list of files
> 
> ```find ./*kallisto21.bam.view.sorted.bam > kallisto21_sorted_bam_list```
> 
> Deduplication using umi_tools
> 
> ```cat kallisto21_sorted_bam_list | parallel -j 12 umi_tools dedup -I {} --output-stats={}.deduplicated -S {}.deduplicated.bam```
> 
> ### 6. Create expression matrix from UMI counts
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







