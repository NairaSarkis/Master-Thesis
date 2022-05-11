# **Naïra Sarkis/Master-Thesis**
## 1. CELSeq2 downstream analysis 

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
> The forward paired reads are used.
> ```demultiplex demux -r -s 1 -e 6 barcodes24.tsv Psam-1_fastp_forward_paired_umiext.fq```
> 
> 
