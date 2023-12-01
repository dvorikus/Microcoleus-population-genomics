# Results from clustering analyses

![clusters](https://github.com/dvorikus/Microcoleus-population-genomics/assets/74075166/d099f208-7c45-4d66-9bab-2e39def53bf6)

Relevant codes for performing fastBAPS (optimized and unoptimized method) and snapclust are in fastBAPS.R

# GTDB-Tk classification

`code()`
conda create -n gtdbtk-2.3.2 conda-forge -c bioconda gtdbtk=2.3.2
conda activate gtdbtk-2.3.2
download-db.sh
gtdbtk classify_wf -x fasta --out_dir gtdbk --genome_dir ./fnas.batch1to7/ --mash_db ./gtdbtk/
