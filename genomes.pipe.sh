#!/bin/bash


# Change the name on line 9
# This script does trimming, assembly, binning, annotation 
# and amphora detection of universal orthologs
# atuhor Petr Dvorak, Palacky University Olomouc
# you can use the script, but there is absolutely nor waranty


# define variable name - sample name, change based on the sample name
READ=A2-D5

# raw data filterind and trimming
java -jar /home/dvorikus/microcoleus.genomes/batch1.2019/trimmomatic/trimmomatic-0.39.jar  PE -threads 24 -phred33  "$READ.1.fastq.gz" "$READ.2.fastq.gz" \
"$READ.output_forward_paired.fq.gz" "$READ.output_forward_unpaired.fq.gz" "$READ.output_reverse_paired.fq.gz" "$READ.output_reverse_unpaired.fq.gz" \
ILLUMINACLIP:/home/dvorikus/microcoleus.genomes/batch1.2019/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

# assembly
spades.py -t 24 \
 --isolate \
 --pe1-1 "$READ.output_forward_paired.fq.gz" \
 --pe1-2 "$READ.output_reverse_paired.fq.gz" \
 --pe1-s "$READ.output_forward_unpaired.fq.gz" \
 --pe2-s "$READ.output_reverse_unpaired.fq.gz" \
 -o spades

# list read files and paths and make read.list
ls -d "$PWD"/*output* > read.list

mkdir maxbin

# copy scaffolds
cp ./spades/scaffolds.fasta ./maxbin

# copy read list
cp read.list ./maxbin

cd ./maxbin 
 
# binning
run_MaxBin.pl -contig scaffolds.fasta -out "$READ.bin" -reads_list read.list -thread 24

# remove scaffold file because they would interfere with the following script
rm scaffolds.fasta

# annotation
for f in *.fasta ; do prokka "$f" --outdir "$f.prokka" ; done

# amphora detection of marker genes
mkdir amphora
cp *.fasta ./amphora
cd ./amphora
for f in *.fasta ; do perl ~/programs/amphora2/AMPHORA2-master/Scripts/MarkerScanner.pl -Bacteria -DNA "$f" ; done

echo done
echo uuuuuuuuuuf
echo computer.tired

