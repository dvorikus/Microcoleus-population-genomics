Here are the files and R scripts used for recombination analysis

# Gubbins

The recent recombination within _Microcoleus_ isolates was estimated with program Gubbins using following command on the alignment generated for the SNP phylogeny (Dataset III)
Alignment available at https://doi.org/10.6084/m9.figshare.24710961.v1

    run_gubbins.py --prefix gubbins_out --first-tree-build rapidnj --tree-builder raxmlng --model GTR alignment.fasta

The outputs are then visualized with Phandango available at http://jameshadfield.github.io/phandango/#/

# Define lineages according to the UPCEL [Kollár et al. 2022](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16218)

![upcel](https://github.com/dvorikus/Microcoleus-population-genomics/assets/74075166/7cc67c74-2841-415f-ba79-0245d17ea2a5)


We have to parse the output from gubbins file recombination_predictions.gff.
From there, we extracted genome fraction that each of the strain shares with other strain belonging to the same lineage (within) and to different lineage (outside) - Supplementary Table 14.
More importantly, we can then extract genome fraction subjected to recombination shared between two lineages as well as genome fraction subjected to recombination shared between strains of one lineage (within) - Supplementary Table 15.

Now we have information on what is the extent of gene flow (homologous recombination) between and within lineage pairs. Using data from Supplementary Table 15, we can calculate genome draction resistant to gene flow between different lineages. This is exactly what we need as a proxy of the divergence probability to define lineages according to the UPCEL.

**Genome fractions shared between lineages (found in Supplementary Table 15), we substract from 100. We also correct the values between the same lineage pairs to 0 (because they are the same lineage).** In Supplementary Table 16 are probabilities of two lineages becoming fully separated species.

According to the UPCEL, these probabilities can range from 0 to 1 (also 0-100%), where values closer to 0 mean that lineages are at early speciation stages and cannot be considered species yet and values closer to 1 mean that lineages have high probabilities of becoming fully separated species.
