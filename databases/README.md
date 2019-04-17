# Default databases

The gene databases used to evaluate the plasmid prediction and binning tool [HyAsP](https://github.com/cchauve/HyAsP).

## MOB-database

Comprises the dereplicated genes of the plasmids listed (by their accession) in `mob_accessions.txt`.
The plasmids and genes had to be at least 500 nt and 100 nt long, respectively.

The database consists of 10591 genes from 230 plasmids of 12 species.

Command:
```
hyasp create mob_database_genes.fasta -a mob_accessions.txt -d -l 500 -m 100
```


## NCBI-database

Comprises the dereplicated genes of the plasmids listed in the plasmid table `plasmids.csv`, 
obtained from [NCBI](https://www.ncbi.nlm.nih.gov/genome/browse#!/plasmids/) on 02 October 2018 (with all possible columns selected via "Choose Columns"), 
but without the plasmids listed in `ncbi_blacklist.txt` or released after 19 December 2015.
The plasmids and genes had to be at least 500 nt and 100 nt long, respectively.
The blacklist contains a few (supposed) *Salmonella enterica* plasmids, 
which have likely been mislabelled as plasmids in NCBI and do not correspond to a *Complete genome*-level assembly.

The database consists of 379917 genes from 5822 plasmids of 1036 species.

Command:
```
hyasp create ncbi_database_genes.fasta -p plasmids.csv -b ncbi_blacklist.txt -d -l 500 -m 100 -t GenBank -r 2015-12-19T00:00:00Z
```

The NCBI-database can also be downloaded from [figshare](https://figshare.com/articles/HyAsP_NCBI_RefSeq_database/8001827).
Note that the database file provided there has to be decompressed before it can be used with **HyAsP**.