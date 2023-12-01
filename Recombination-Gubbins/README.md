Here are the files and R scripts used for recombination analysis

# Gubbins

The recent recombination within _Microcoleus_ isolates was estimated with program Gubbins using following command on the alignment generated for the SNP phylogeny (Dataset III)
Alignment available at https://doi.org/10.6084/m9.figshare.24710961.v1

    run_gubbins.py --prefix gubbins_out --first-tree-build rapidnj --tree-builder raxmlng --model GTR alignment.fasta

The outputs are then visualized with Phandango available at http://jameshadfield.github.io/phandango/#/
