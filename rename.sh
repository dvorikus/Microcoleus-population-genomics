#!/bin/bash

## parse results from Snippy - consensus.subs.fa
## all samples are concatenated into one huge alignment

# remove plasmid names in file
for f in *.fa
do
	sed -i 's/>CP00361[5-9].1//g' "$f"
done

# remove blank lines within the file
for f in *.fa
do
	sed -i '/^\s*$/d' "$f"
done


for f in *.fa
do
    name="$f"
    sed -i s/CP003614.1/"$name"/ "$f"
done

for f in *.fa
do
	sed -i 's/.output_forward_paired.fq.gz.consensus.subs.fa//g' "$f"
done


cat *.fa > all.samples.fa

echo done
