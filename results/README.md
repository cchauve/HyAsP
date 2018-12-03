# Overview

Here, we provide the results of an comparative evaluation of **HyAsP** with plasmidSPAdes and MOB-recon 
on a collection of real plasmids originally compiled for benchmarking [MOB-suite](https://dx.doi.org/10.1099/mgen.0.000206).
The benchmarking data set consists of 133 bacterial samples comprising the same number of closed genomes and 377 plasmids (`samples.csv`). 
In order to simulate the use of plasmids identification tools on a newly WGS dataset, using only knowledge from previously existing plasmid data, 
the data set was split in half based on the release date of the corresponding read data.
The 66 samples released on 19 December 2015 or later formed the *test samples* (listed in `test_ids.txt`).
We assessed the predictions of the three tools on each test sample in terms of their precision, recall and F1 score by mapping the predictions against the expected plasmids.
 
The other 67 samples of the benchmarking data set were used as the first set of *database samples* and we refer to the resulting database as the **MOB-database**.
The results of this assessment are stored in `analysis_mob_filtered/` and all scores are compiled in `analysis_mob_filtered/scoring_results.csv`. 
For each sample, `analysis_mob_filtered/` also has a correspondingly named subdirectory which contains a results file
for the comparison between the predicted and expected plasmids per tool.
For example, the results file for **HyAsP** on sample 1 is stored in `analysis_mob_filtered/sample_1/eval/greedy/greedy_eval.csv`.
The results for all test samples were evaluated and visualised in `evaluate_mob_filtered.ipynb`.   
In addition, we examined the occurrence of misassembly events in the plasmid sequences produced by **HyAsP** through QUAST analyses.
For each sample, `quast_results_mob_filtered/` contains a subdirectory providing the overall QUAST report (`report.tsv`) and the misassembly report (`contigs_reports/misassemblies_report.tsv`).
In both reports, the columns *greedy* and *greedy_contigs* refer to **HyAsP** when using the plasmid sequences respectively the underlying contig collections in the QUAST analysis, wihout accounting for the order in which they were chained.
The QUAST results were summarised in `quast_analysis_mob_filtered.ipynb`. 
A comparison with plasmidSpades and MOBrecon is provided, however both methods, unlike **HyAsP** provide only collections of contigs, that can then be freely positioned and orientated by QUAST and, thus, show rarely any misassembly event (which, in this case, originate from the intial read assembly). 

Similarly, we assessed the tools on a second set of database samples based all plasmids available from NCBI that were released before 19 December 2015 (**NCBI-database**). 
The analysis files and results (in `analysis_ncbi_filtered` etc.) have the same structure as for the analysis using the MOB-database.
