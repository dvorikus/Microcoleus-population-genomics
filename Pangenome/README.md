# Coinfinder

![coinfinder](https://github.com/dvorikus/Microcoleus-population-genomics/assets/74075166/3210fa44-fcdb-46ad-a3b0-c551c3d086f4)

Coinfinder was used to determine the gene cooccurrence and avoidance within _Microcoleus_ flexible genome. As an input, we used gene presence/absence matrix computed from Roary. List of statistically significant associating and dissasociating genes is in Supplementary Table 7.

    coinfinder -I gene_presence_absence.csv -p core_gene_alignment.aln.treefile -o coinfinder.associate --associate -m
    coinfinder -I gene_presence_absence.csv -p core_gene_alignment.aln.treefile -o coinfinder.dissassociate --dissociate -m


# HGTector2

We investigated signature of horizontal gene transfer (HGT) using HGTector diamond all-against-all similarity search on dataset III. As an input, we used multifasta files of AA sequences (201 genomes, without the outgroup M2_D5).

First we install the program if we do not have it

    conda create -n hgtector -c conda-forge python=3 pyyaml pandas matplotlib scikit-learn bioconda::diamond
    conda activate hgtector
    pip install git+https://github.com/qiyunlab/HGTector.git

Then, we downloaded all protein sequences of NCBI RefSeq genomes of bacteria (322,197), archaea (1866), fungi (551), and protozoa (96), keeping 1 genome per species.

    hgtector database -o ./hgtector/db/ --default

Now, we can perform the search for potentially horizontally transfered genes

    hgtector search -i /media/second/microcoleus/faas/ -o ./search/ -t ./db/taxdump/ -d ./db/diamond/db.dmnd -m diamond -p 16

Finally, we input results from previous search and manually parse and inspect files

    hgtector analyze -i ./search/ -o ./analyze/ -t ./db/taxdump/ --donor-name 

_Note_ All the files necessary for making Figure 3 can be found here.
