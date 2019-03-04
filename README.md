# HyAsP

**HyAsP** (Hybrid Assember for Plasmids) is a tool for the extracting plasmids from WGS assemblies in an automatic way.
It combines ideas from both reference-based and depth-based methods to identify plasmids in a greedy algorithm, 
using information on the occurrences of known plasmid genes and considering characteristics of the contigs such as read depth and GC content.



### Overview

Directory `HyAsp/` contains the source code of **HyAsP**, which can be installed as a package through `setup.py` (see below).
Directory `databases/` provides exemplary files that can be used to construct a gene database for **HyAsP**, 
while `results/` contains the results of a comparison of **HyAsP** with plasmidSPAdes and MOB-recon. 



## Requirements

**HyAsP** was developed and tested with the following software dependencies: 
  - [Python](https://www.python.org/downloads/) (`python`, version 3.5.4; *packages*: Bio, math, numpy, os, pandas, random, subprocess, sys)
  - [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (`makeblastdb` `tblastn` and `blastn`; version 2.6.0)
  - standard UNIX tools (`mkdir`, `rm`, `cat`)
  
BLAST+ is only required for the `create` command.
  
In order to use **HyAsP** as part of the provided pipeline starting from FASTQ reads, 
the following requirements have to be satisfied in addition:
  - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (`fastqc`, version 0.11.5)
  - [sickle](https://github.com/najoshi/sickle) (`sickle`, version 1.33)
  - [cutadapt](https://cutadapt.readthedocs.io/en/stable/) (`cutadapt`, version 1.16)
  - [Trim Galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) (`trim_galore`, version 0.4.5_dev)
  - [Unicycler](https://github.com/rrwick/Unicycler) (`unicycler-runner.py`, version 0.4.5)
    - [SPAdes](http://cab.spbu.ru/software/spades/) (`spades.py`, version 3.12.0)
    - [Racon](https://github.com/isovic/racon) (`racon`, version 1.3.0)
    - [Pilon](https://github.com/broadinstitute/pilon/wiki) (`pilon-1.22.jar`, version 1.22)
    - [SAMtools](http://www.htslib.org/) (`samtools`, version 1.5)
    - [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (`bowtie2` and `bowtie2-build`; version 2.3.3.1)
    - [Java](https://www.java.com/en/download/manual.jsp) (`java`, version 1.8.0_121)
 
Python and the other tools have to be in the `PATH` or specified through their path options.



## Installation

Get the source code from GitHub and, optionally, install **HyAsP** as a package:
```
git clone https://github.com/romueller/hyasp.git
cd hyasp
python setup.py sdist
pip install dist/HyAsP-1.0.0.tar.gz
```
Installing **HyAsP** as a package makes `hyasp.py` and `fastq_to_plasmids.py` available as `hyasp` and `hyasp_pipeline`, respectively.
We strongly recommend installing **HyAsP** in a virtual environment which uses Python 3.

Subsequently, the proposed default gene database could be built:
```
hyasp create databases/default_genes_db.fasta -d -p databases/plasmids.csv -b databases/default_blacklist.txt
```
The plasmid table (`databases/plasmids.csv`) was downloaded on 21 November 2018 from 
[NCBI](https://www.ncbi.nlm.nih.gov/genome/browse#!/plasmids/) (with all possible columns selected via "Choose Columns")
and reformatted by changing the separator to tab and removing quotes surrounding values.
The blacklist (`databases/default_blacklist.csv`) contains a few (supposed) *Salmonella enterica* plasmids,
which have turned out to be problematic (i.e. too chromosome-like) in our analyses.
However, the gene database does not have to be created using the `create` command of **HyAsP** or with above parameters.



## Usage

**HyAsP** provides the functions to find plasmids and create necessary inputs through different commands.
Below example show simple uses of the commands.
See section *Parameters* for lists of options to change the behaviour of each command.

 
 
### 1) Create a gene database from a collection of plasmids

Command: `hyasp create`

Requires a list of plasmids with annotated genes or a plasmid table obtained from [NCBI](https://www.ncbi.nlm.nih.gov/genome/browse#!/plasmids/) (with all columns).

**Example:**
```
hyasp create genes.fasta -a accessions.txt
```
Here, the gene database (`genes.fasta`) is created by downloading the GenBank files of the plasmids given by accession number in `accessions.txt` (one per line)
and extracting the genes from them.

**Additional options:** 
```
--from_accession, -a        Path to file containing one (plasmid) accession number per line OR list of accession numbers (see --from_command_line).
                            (default: (empty string), i.e. not used)
--from_genbank, -g          Path to file containing path to a GenBank file per line OR list of paths (see --from_command_line).
                            (default: (empty string), i.e. not used)
--from_plasmid_table, -p    Path to plasmid table downloaded from NCBI.
                            (default: (empty string), i.e. not used)
--keep_plasmids, -k         Stores the plasmids underlying the gene database in FASTA format if a file is specified.
                            (default: (empty string), i.e. deactivated)
--dereplicate, -d           Removes duplicate genes from database if activated.
                            (default: False)
--from_command_line, -c     Instead of a file containing the accession numbers (file paths), the options -a (-g) expect
                            a comma-separated list of accession numbers (file paths). Cannot be combined from -p.
                            (default: False)
--extend, -e                Genes (and plasmids) are added to an existing database instead of overwriting it.
                            (default: False)
--released_before, -r       Consider only plasmids released before the specified date. Can only be combined with -p.
                            Date format: YYYY-MM-DDTHH:MM:SSZ, e.g. 2005-07-31T00:00:00Z.
                            (default: (empty string), i.e. deactivated)
--type, -t                  Build the databases from the RefSeq accession numbers (RefSeq), GenBank accession numbers (GenBank) or both (both).
                            Affects only option -p.
                            (default: both)
--blacklist, -b             Comma-separated list of accession numbers of plasmids not to be included in the databases.
                            Cannot be combined with -e.
                            (default: (empty string), i.e. deactivated)
--min_length, -l            Minimum length of plasmids to be considered for the database.
                            (default: 0)
--max_length, -L            Maximum length of plasmids to be considered for the database.
                            (default: infinity)
```
The database is created from either accession numbers or (already downloaded) GenBank files or the NCBI plasmid table, i.e. the options `-a`, `-g` and `-p` cannot be combined.



### 2) Map a collection of genes to the contigs of an assembly

Command: `hyasp map`

Requires a gene database and a collection of contigs.

**Example:**
```
hyasp map genes.fasta gcm.csv -g assembly.gfa  
```
Here, the gene-contig mapping `gcm.csv` is determined by mapping the genes in `genes.fasta` to the contigs of the assembly
in `assembly.gfa`. 

**Additional options:**
```
--from_fasta, -f    Path to the file containing the contigs (in FASTA format) to which the genes should be matched.
                    (default: (empty string), i.e. not used)
--from_gfa, -g      Path to the file containing the contigs (as part of an assembly graph in GFA format) to which the genes should be matched.
                    (default: (empty string), i.e. not used)
--clean, -c         Remove temporary files after the mapping has been created.
                    (default: False)
--makeblastdb       Path to the makeblastdb executable.
                    (default: makeblastdb)
--blastn            Path to the blastn executable.
                    (default: blastn)
```
The contigs are read from either a FASTA or a GFA file, i.e. either `-f` or `-g` have to be used.



### 3) Filter a gene-contig mapping

Command: `hyasp filter`

Requires a gene database and a gene-contig mapping (created with `map`).

**Example:**
```
hyasp filter genes.fasta gcm.csv filtered_gcm.csv 
```
Here, the gene-contig mapping `gcm.csv` (based on the genes in `genes.fasta`) will be filtered using default thresholds 
(for identity and length) and the remaining mapping is stored in `filtered_gcm.csv`. 

**Additional options:**
```
--identity_threshold, -i    Minimum identity of hits retained in the mapping.
                            (default: 0.95)
--length_threshold, -l      Minimum fraction of query (gene) that has be matched to keep a hit.
                            (default: 0.95) 
--find_fragmented, -f       Search for fragmented hits, i.e. several short high-identity hits that together satisfy the length threshold.
                            (default: False)
```



### 4) Find plasmids in an assembly graph

Command: `hyasp find`

Requires an assembly graph (in [GFA v1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) format), 
a gene database (see `create`) and a (filtered) gene-contig mapping (`map`, `filter`).

**Example:**
```
hyasp find assembly.gfa genes.fasta gcm.csv output_dir
```
Here, the plasmids are predicted for the assembly provided in `assembly.gfa`, based on the gene-contig mapping `gcm.csv`
and the gene database provided in `genes.fasta`.

**Additional options:**
```
--min_gene_density, -g          Minimum gene density of a putative plasmid. Plasmids with a lower gene density are marked as questionable.  
                                (default: 0.3) 
--min_seed_gene_density, -k     Minimum gene density necessary for a contig to be considered as a seed.     
                                (default: 1.5 * min_gene_density) 
--min_length, -l                Minimum length of a putative plasmid. Shorter plasmids are marked as questionable.
                                (default: 1500) 
--max_length, -L                Maximum length of a putative plasmid. Gene-containing contigs longer than max_length are not used as seeds. 
                                A contig is excluded from list of potential extensions, if the combined length of the contig 
                                and the plasmid is larger than max_length.
                                (default: 1750000) 
--min_read_depth, -r            Minimum read depth of a contig to be able to participate in a plasmid.
                                (default: 0.75 * min_plasmid_read_depth) 
--min_plasmid_read_depth, -d    Minimum average read depth of a putative plasmid. Plasmids with a lower average read depth are marked as questionable.
                                (default: 0.4 * (median read depth of input assembly graph)) 
--max_gc_diff, -G               Maximum difference in GC content between a plasmid and a potentially added contig.
                                (default: 0.15) 
--max_intermediate_contigs, -c  Maximum number of gene-free contigs between two gene-containing contigs in a plasmid. 
                                A contig is excluded from the list of potential extensions, if its addition would violate this threshold. 
                                (default: 2) 
--max_intermediate_nt, -n       Maximum total length of any consecutive sequence of gene-free contigs in a plasmid. 
                                A contig is excluded from the list of potential extensions, if its addition would violate this threshold. 
                                (default: 2000) 
--max_score, -s                 Maximum score of a potential extension. Possible extensions with a higher score are discarded. 
                                (default: infinity) 
--score_weights, -w             Weights of the different components of the function used to score extensions. 
                                Comma-separated list of entries of the form <name>=<value>.
                                The weight of a component is determined automatically if the question mark (?) is used as <value>. 
                                (default: depth_diff=1,gene_density=1,gc_diff=1) 
--keep_subplasmids, -q          Do not mark plasmids whose underlying set of contigs is contained by others as questionable. 
                                If several plasmids have the same underlying set of contigs, one of them will remain 
                                in the collection of putative plasmids. 
                                (default: False) 
--overlap_ends, -o              Minimum overlap between the two ends of a plasmid in order to mark it as circular. 
                                (default: infinity) 
--binning, -b                   Factor determining how many standard deviations the read depth and GC content of plasmids are allowed to differ
                                from the 'centre' of their bin. Binning is activated by setting the parameter to any value different from NaN. 
                                (default: NaN, i.e. deactivated) 
--fanout, -f                    Maximum number of predecessors / successors of any contig in a 'plasmid' (or rather contig collection). 
                                Setting this parameter to any value > 1 leads to non-linear / branching contig chains. 
                                Changes which files are generated as output. Cannot be used together with probabilistic.
                                (default: 1) 
--probabilistic, -p             Flag changing the behaviour of the extension step to a probabilistic choice. 
                                The probability of an extension is the share of involved contig of the total read depth of all extensions. 
                                Cannot be used together with fanout. 
                                (default: False) 
--use_node_based, -N            Flag changing the behaviour of the extension step to a node-based loop avoidance (instead of link-based).           
                                (default: False) 
--use_median, -u                Flag activating the use of median (instead of mean) in order to compute the average read depth of a plasmid. 
                                (default: False) 
--verbose, -v                   Flag activating detailed logging (of the extension procedure). 
                                (default: False) 
```



## Outputs

**Contig chains**  
Lists contigs and their orientation as they appear in the linear contig chain of all plasmids (both putative and questionable).
Only for `fanout = 1`.    
*File name:* `contig_chains.csv`    
*Format:* `<plasmid id>;<comma-separated list of contigs with orientation>`    
*Example:* plasmid_0;23+,25-,10+    

**Plasmids**    
Stores the plasmid sequences (concatenations of the (orientated) contigs) in FASTA format. 
The plasmid identifier is also used as the identifier of the FASTA entry. 
The additional information on each plasmid are provided in the deflines (e.g. seed_contig and gene_density.
Only for `fanout = 1`.        
*File names:* `putative_plasmids.fasta`, `questionable_plasmids.fasta`    
*Format:* FASTA with additional information in defline (as tab-separated list of `<property>=<value>` pairs)

**Contig collections**    
Lists name and orientation of all contigs for each putative resp. questionable plasmid.
Only for `fanout > 1`.    
*File names:* `putative_plasmid_contigs_list.csv`, `questionable_plasmid_contigs_list.csv`     
*Format:* `<plasmid id>;<comma-separated list of contigs with orientation>`    
*Example:* `plasmid_0;23+,25-,10+`

**Contigs**    
Stores the contigs underlying plasmids in FASTA format. 
If a contig is used in negative orientation, the reverse complement of its sequence is stored in the output file.
The identifier of a FASTA entry consists of the contig name and the contig identifier (separated by the `|`-symbol).    
*File names:* `putative_plasmid_contigs.fasta`, `questionable_plasmid_contigs.fasta`        
*Format:* FASTA    

**Tagged assembly graph**    
Stores a copy of the input assembly graph and adds colour and label information to contigs used in putative and questionable plasmids.
Contigs occurring in at least one putative plasmids are blue, while those occurring one or more questionable plasmids 
(but no putative plasmid) are light blue.
Each contig is labelled with the identifiers of the plasmids it occurs in and seed contigs also contain a `*` in their label.    
*File name:* `tagged_assembly.gfa`    
*Format:* GFA with additional (optional) tags for contigs


**Plasmid bins**    
Lists the plasmid identifiers (of putative plasmids resp. all plasmids) grouped into the different bins.
Only for `binning != NaN`.    
*File names:* `plasmid_bins_putative.csv`, `plasmid_bins_all.csv`    
*Format:* `<comma-separated list of plasmid identifiers>`    
*Example:* `plasmid_0,plasmid_10,plasmid_3`    



## Pipeline from FASTQ reads to plasmids

The pipeline takes FASTQ reads and a collection of (plasmid) genes as input. 
The reads are preprocessed (sickle, Trim Galore) and assembled (Unicycler).
The preprocessing and an analysis of the read data (FastQC) are optional.
Subsequently, the genes are mapped to the assembly contigs using BLAST (blastn / megablast) 
and plasmids are predicted in the assembly graph using **HyAsP**.



### Usage

The pipeline can be used as a stand-alone script or imported as a module.

The simplest usage only requires the output directory, a gene database and (unpaired) short FASTQ reads.

```
hyasp_pipeline output_dir genes.fasta -s reads.fastq
```


Paired short-read data can be used by using the `-1` and `-2` options, e.g.
```
hyasp_pipeline output_dir genes.fasta -1 first_reads.fastq -2 second_reads.fastq 
```

Long reads can be added (to unpaired and / or paired short-read data) through the `-l` option, e.g.
```
hyasp_pipeline output_dir genes.fasta -s short_reads.fastq -l long_reads.fastq
```

Unicycler's assembly mode can be changed from `normal` via the `-u` option and another gene database can be specified (`-D`), e.g.
```
hyasp_pipeline output_dir genes.fasta -s reads.fastq -u conservative
```

The gene-contig mapping obtained from BLAST is filtered before the greedy algorithm is used.
The filtering can be influenced by changing the length and / or identity threshold, e.g.
```
hyasp_pipeline output_dir genes.fasta -s reads.fastq --identity_threshold 0.9 --length_threshold 0.92
```
The default value for both is 0.95.

In addition, the options of the greedy algorithm can be given to `hyasp_pipeline`. 
`hyasp_pipeline -h` provides a list of the possible options.


### Outputs

The final output of the pipeline corresponds to the output of **HyAsP** described above and is stored in the `plasmids/` subdirectory of the output directory.

The results of the preprocessing and assembly step are found in the `data/` and `assembly/` subdirectory of the output directory, respectively.

