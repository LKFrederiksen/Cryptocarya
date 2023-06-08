# Supplementary information for the article manuscript ["Long-distance dispersal shaped the pantropical distribution of the geographically structured Cryptocarya group"]() from the MSc thesis of Laura Kragh Frederiksen.

##**DATA**

**``Data_S1.zip``**

**Data S1: Alignments**
342 fasta files generated with MAFFT and trimmed with Trimal, Optrimal, CIalign and Taper. Each file represents one gene retrieved from the ['mega353' target file](https://github.com/chrisjackson-pellicle/NewTargets) (McLay et al., 2021). The faster headers correspond to each individual in short notation. See attached [text file]() for a conversion of the short notation to the species name.

**``Data_S2.zip``**

**Data S2: Gene trees**
342 ML (maximum likelihood) gene trees generated with IQtree in Newick format. Each gene tree corresponds to one gene. The faster headers correspond to each individual in short notation. See attached [text file]() for a conversion of the short notation to the species name.

##**FIGURES**

**``Figure_S1.pdf``**

**Figure S1: Dispersal rate scalers for BioGeoBEARS** 
Adapted from [Toussaint et al. (2017)](https://onlinelibrary.wiley.com/doi/10.1111/jbi.12977). Dispersal rate scaling matrices used in the [BioGeoBEARS](https://doi.org/10.21425/F5FBG19694) (Matzke, 2013) analyses from the oldest time slice to the newest (TS1 to TS3). A matrix of dispersal rate scalers relative to potential geographic barriers is given for each of the three time slices: ~80 to 40 Ma, 40 to 20 Ma and 20 Ma to present. Four symbols are used to indicate the kind of barrier between the selected ranges. Regions are abbreviated as follows: A = Andean-Argentinian, B = Neotropical, C = Southern Africa, D = African, E = Madagascan, F = Northern Australia, G = Malesian, H = Indian-Indochinese, I = Neozealandic-Patagonian and J = Eurasiatic (see [Carta et al., 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9298788/)).

**``Figure_S2.pdf``** 

**Figure S2: Species tree with quartet scores** 
[ASTRAL](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5998893/) (Zhang et al., 2018) species tree including ingroup and outgroup for all samples (also species with multiple samples) with quartet scores (QS). The quartet scores are included as pies on each node. Red = main topology, turquoise = alternative topology, grey = second alternative topology.

**``Figure_S3.pdf``**

**Figure S3: Biogeograchic ancestral range estimation using BioGeoBEARS including range probability scores**
Ancestral range estimation for the Cryptocarya group using BioGeoBEARS based on the DEC+TS+_j_ model, with ancestral ranges, and [dispersal rate scalers]() between regions for each time slice (TS1, TS2, TS3). Roman numerals on the nodes represent crown nodes of important colonisation events, discussed in the text and shown in [Table 4](). Pies on the nodes and branches show the probality of each range being the 'true' ancestral range. These ranges can comprise of one to three combined regions, as indicated by single letters or a combination of these (e.g., FG = Northern Australia and Malesia). Regions are abbreviated as follows: A = Andean-Argentinian, B = Neotropical, C = Southern Africa, D = African, E = Madagascan, F = Northern Australia, G = Malesian, H = Indian-Indochinese, I = Neozealandic-Patagonian and J = Eurasiatic (see [Carta et al., 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9298788/). Three modes of range inheritance (range expansion, local extinction, and narrow vicariance) are indicated as symbols on branches in the phylogeny. 

**``Figure_S4.pdf``**

**Figure S4: Dispersal sinks and sources - BioGeoBEARS BSM**
Number of dispersal events (range expansion and jump dispersal) for each time slice (TS1, TS2, TS3; see Figure 3 in [thesis]()) from the DEC+TS+_j_ model (see [Table 3]()). Counts of events were averaged across the 100 biogeographic stochastic mappings (BSMs) with standard deviations in parentheses. Rows represent source ranges; columns represent dispersal sinks. Darker shades indicate a higher frequency of dispersal events. The sum and percent of events in each row and column are given on the margins. Regions are abbreviated as follows: A = Andean-Argentinian, B = Neotropical, C = Southern Africa, D = African, E = Madagascan, F = Northern Australia, G = Malesian, H = Indian-Indochinese, I = Neozealandic-Patagonian and J = Eurasiatic.

##**TABLES**

**``Table_S1.xlsx``**

**Table S1: Gene tree concordance to 'PAFTOL' outgroup topology**
Nuclear [Angiosperms353](https://www.ncbi.nlm.nih.gov/pubmed/30535394) (Johnson et al., 2019) genes used to generate the [PAFTOL](https://www.ncbi.nlm.nih.gov/pubmed/33983440) (Baker et al., 2022) outgroup topology. We used [SortaDate](https://www.ncbi.nlm.nih.gov/pubmed/29772020) (Smith et al., 2018) to identify the supercontig genes that were most concordant to the PAFTOL outgroup topology based on exons (Baker et al., 2022). We then tested how many exons we were able to substitute with supercontigs without losing the desired topology used in [Ramirez-Barahona (2020)](https://www.nature.com/articles/s41559-020-1241-3). This was 10%. The genes in the supercontig table are the 10% most concordant genes, whereas the genes in the exon table are the remaining 90 %.

**``Table_S2.xlsx``**

**Table S2: HybPiper recovery stats**
[HybPiper](https://bioone.org/journals/applications-in-plant-sciences/volume-4/issue-7/apps.1600016/HybPiper--Extracting-Coding-Sequence-and-Introns-for-Phylogenetics-from/10.3732/apps.1600016.full) target enrichment and gene recovery efficiency summary for the 115 ingroup samples (*Syndiclis chinensis* included after the ['coverage']() step and 26 outgroup samples (*Peumus boldus* and *Chimonanthus salicifolia* included after ['coverage']() step). The columns in the table correspond to: Sample name , Number of reads,	Number of reads on target, Percentage of reads on target,	Number of genes with reads,	Number of genes with contigs,	Number of genes with sequences,	Number of genes with sequences > 25% of the target length, Number of genes with sequences > 50% of the target length, Number of genes with sequences > 75% of the target length, Number of genes with sequences < 150% of the target length, Number of genes with paralog warnings. See [HybPiper wiki](https://github.com/mossmatters/HybPiper/wiki#hybpiper-stats) for an explanation of what the different outputs mean. See attached [text file]() for a conversion of the short notation to the species name.

**``Table_S3.xlsx``**

**Table S3: Paralog analyses from HybPiper**
[HybPiper](https://bioone.org/journals/applications-in-plant-sciences/volume-4/issue-7/apps.1600016/HybPiper--Extracting-Coding-Sequence-and-Introns-for-Phylogenetics-from/10.3732/apps.1600016.full) generated output showing potential paralogs for each gene and number of gene duplications for each sample. 0 indicates no genes recovered for the specific sample, 1 indicates gene recovered for the sample, >1 indicated that the gene potentially is paralogous. 

**``Table_S4.xlsx``**

**Table S4: BioGeoBEARS dispersal sinks and sources**
Number of dispersal (range expansion, e; jump dispersal, _j_) events for the three time slices (see [Figure 3](); 80 Ma, 40 Ma, 20 Ma and 0 Ma) from the DEC+TS+j model (see [Table 3]()). Counts of events were averaged across the 100 biogeographic stochastic mappings (BSMs) with standard deviations in parentheses. Rows represent source ranges; columns represent dispersal sinks. Darker shades indicate a higher frequency of dispersal events. The sum and percent of events in each row and column are given on the margins. Regions are abbreviated as follows: A = Andean-Argentinian, B = Neotropical, C = Southern Africa, D = African, E = Madagascan, F = Northern Australia, G = Malesian, H = Indian-Indochinese, I = Neozealandic-Patagonian and J = Eurasiatic.

##**TEXT**

**``Text_S1.pdf``**

**Text S1: Fossil selection**
Justification for the eight fossils included as calibration points in our molecular dating with [BEAST](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6007674/) based on [Ramirez-Barahona et al. (2020)](https://www.nature.com/articles/s41559-020-1241-3) and [Li et al. (2020)](https://www.ncbi.nlm.nih.gov/pubmed/32619568).
