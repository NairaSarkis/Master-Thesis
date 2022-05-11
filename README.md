# **Naïra Sarkis/Master-Thesis**
## 1. P.sambesii whole genome annotation
### 1.1 Obtain RNA-seq data from NCBI and prepare read files
> Use NCBI-SRA-toolkit to download paired-end RNA-seq data for P.sambesii that has been generated before (###Ref)
> 
> ```fasterq-dump -S SRR8243961```
> 
> Headers of read files are modified using sed. Dots in headers are exchanged with underspace and everything starting from the first blank space is deleted.
> 
### 1.2 Mask the genome using repeatmodeler and repeatmasker
> ```BuildDatabase -name ES601_gene_DB -engine ncbi psambesii_genome.fasta RepeatModeler -engine ncbi -pa 16 -database ES601_gene_DB```
>
> ```RepeatClassifier -consensi ES601_gene_DB-families.fa```
>
> ```RepeatMasker -pa 16 -e ncbi -lib ES601_gene_DB-families.fa.classified psambesii_genome.fasta```
### 1.3 Index genome and align reads using gmap/gsnap
> ```gmap_build -D /scratch/nsarkisk/Psam_annotation -d genome_index psambesii_genome.fasta.masked```
> 
> ```gsnap -D /scratch/nsarkisk/Psam_annotation -d genome_index -A sam -o /scratch/nsarkisk/Psam_annotation/psambesii-gsnap.sam SRR8243961_1.sednew.fastq SRR8243961_2.sednew.fastq```
### 1.4 Gene predictions using the braker2 pipeline
> ```braker.pl --species=PlectusSambesii --softmasking --AUGUSTUS_CONFIG_PATH=/scratch/nsarkisk/Psam_annotation/augustus-config/ --genome=psambesii_genome.fasta.masked --bam=psambesii-gsnap.bam.sorted```

## 2. CELSeq2 downstream analysis 

### Tools implemented in this analysis:

- Table of versions
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
> Complexity describes the percentage of base that is different from its next base. The complexity filter default value is 30, which means 30 % complexity is required, to keep the reads. (SOURCE?)
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
> --------------------------------------------------TABLE------------------------------
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
> Indexing the genome using kallisto
> ```kallisto index -i P_sambesii21.index -k 21 psambesii_genome.fasta.masked```
