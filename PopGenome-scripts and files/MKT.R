#Perform McDonald-Kreitman test for positive selection with PopGenome

library(PopGenome)
microcoleus_vcf <- readVCF("./vcf/all.biallelic.version2.02.vcf.gz", tid = "CP003614.1", frompos = 1, topos = 7479014, include.unknown = TRUE, numcols = 2000, gffpath = "./gff/reference_ext_tab.gff")
#gff file has to have the same chromosome name
populations_microcoleus_vcf <- read.csv("populations_microcoleus_vcf.txt", sep="")
populations_vcf <- split(populations_microcoleus_vcf$ind, populations_microcoleus_vcf$pop)
microcoleus_vcf <- set.populations(microcoleus_vcf, populations_vcf, diploid = F)

#Here we need a full genome of our reference genome Oscillatoria nigro-viridis PCC 7112 obtained from the GenBank
microcoleus_vcf <- set.synnonsyn(microcoleus_vcf, ref.chr = "./GCF_000317475.1_ASM31747v1_genomic.fna" )

#You can just remove the sequence bellow annotation (dont forget to hit Enter and go to the new row)
#Be wary that sometimes get_gff_info extracts incomplete list of genes present in annotation. If that occurrs, usually it's due to the formatting of the gff file
genePos <- get_gff_info(gff.file = "gff/Reference_ext_tab.gff", chr = "CP003614.1", feature = "CDS")
genes <- splitting.data(microcoleus_vcf, positions = genePos, type = 2)
split <- MKT(genes, do.fisher.test = TRUE)

#to name the regions by genes names
gene_name <- array()
for (j in 1:length(split@region.names)){
  gene_name[j] <- head(get_gff_info(split, position=j, chr="CP003614.1", gff.file="./gff/reference_ext_tab.gff"))
}
#if this doesnt work, use shell after the next step and parse the files to get proper gene names for each pairwise MK test

for(i in 1:nrow(split@MKT)){
     write.csv(split@MKT[i,], file=paste0(i, ".csv"))
}

#After this step, we can also use shell to extract only statistically significant MK results that indicate positive selection
#We can also manually check which genes are under positive or negative selection, by going through the output files from the previous command
#We will have one csv file per gene named numerically from 1, which correspond to the 1st gene in the annotation file.