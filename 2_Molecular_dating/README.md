# Dating the phylogeny with BEAST

[BEAST](https://beast.community) [(Drummond & Rambaut, 2007)](https://www.ncbi.nlm.nih.gov/pubmed/17996036
https://bmcecolevol.biomedcentral.com/counter/pdf/10.1186/1471-2148-7-214.pdf) based on Bayesian inference was used to date the phylogenetic tree. To do this we used the three best nuclear genes, nine fossils (primary calibrations) and one secondary calibration on the root (see Table ... in [MSc thesis]()). Moreover, we fixed the tree topology to that of the ASTRAL tree based on 34. nuclear genes. 

We selected the fossils based on [Ramirez-Barahona et al., 2022](https://www.nature.com/articles/s41559-020-1241-3) and the fossils included in [Li et al. (2020)](https://www.ncbi.nlm.nih.gov/pubmed/32619568), except the fossil [_Cryptocaryoxylon gippslandicum_](https://doi.org/10.1080/03115518608619157).

Look at part 3 of my [workflow.py]() file based on [gwf]() to see the steps and code I used to undertake parts of the molecular dating analysis.

## Selecting genes

To select genes for dating I used [SortaData](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0197433&type=printable) to get a list of the 'best' genes, thus sorted after the genes with least root to tip variance and best bipartition support compared to the species tree. See [workflow.py](). This is done in order to select the genes, whose gene trees have the least topological conflict with the species tree. Here the three genes (5594, 5620 and 6139) were selected.

## Alignments
Based on the output of SortaDate I selected the three 'best' genes that contained 'supercontigs' (exons including the flanking intron and spacer-regions) in both the ingroup and outgroup. 
For each gene in the alignments, if a sample dit not have a sequence, as no reads had been recovered for the specific gene, then gaps were added as required for BEAST to run. This was done using [script]() based on Biopython v1.7.5, an adapted version of a script [/1.1_add_empty_seqs_for_sp.py](2_phylo_dating/1.1_add_empty_seqs_for_sp.py)created by Maya Schrödl.

For the species where more than one sample was available, only one individual was chosen to comply with BioGeoBEARS. The individual with the alignment that had the least gaps was selected.
These genes were manually cleaned using the program [](). 

The cleaned alignments can be retrieved from the folder []().

The tips of the samples that were in excess were dropped from the ASTRAL species tree, so it could be used as starting tree and to fix the topology. The starting tree can be retreived here: [tree]().

## Preparing the xml file in BEAUti

To run the molecular dating analysis in BEAST a xml file is required. This file can be generated in the BEAST GUi [BEAUti](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3408070/pdf/mss075.pdf). I used version 1.10.4. 

The best-performing evolutionary model for each gene was determined using the Akaike Information Criterion (AIC; Akaike, 1974) as implemented in IQ-TREE ModelFinder [(Kalyaanamoorthy et al., 2017)](https://www.nature.com/articles/nmeth.4285.pdf). The GTR+Γ+I with 4 Γ categories was selected for gene 5620 and 5594, whereas HKY+ Γ with 4 Γ categories was selected for gene 6139. 

An uncorrelated relaxed clock model assuming an uncorrelated lognormal rate prior and a birth-death speciation prior was selected. The ASTRAL tree was selected as user-specified starting tree and used as constraint by deleting the commands: subtreeSlide; narrowExchange; wideExchange; wilsonBalding in BEAUti v.1.10.4, because a topology inferred from 300+ genes is more reliable than one inferred from three.

### Steps in BEAUti for each tab

`Partitions`
<img width="1067" alt="image" src="https://github.com/LKFrederiksen/Cryptocarya/assets/112875495/50ebf68c-7f1a-4a03-bdf0-892dc8255354">

`Taxa`
<img width="1061" alt="image" src="https://github.com/LKFrederiksen/Cryptocarya/assets/112875495/a719ed70-cf59-4064-bc44-51054d4d8dd7">

`Sites`
<img width="1066" alt="image" src="https://github.com/LKFrederiksen/Cryptocarya/assets/112875495/71ef6058-d718-4b19-be37-40be23a13cf2">

`Clocks`
<img width="1066" alt="image" src="https://github.com/LKFrederiksen/Cryptocarya/assets/112875495/46e3074e-0034-4ec7-ae0f-437e517367ed">

`Trees`
<img width="1066" alt="image" src="https://github.com/LKFrederiksen/Cryptocarya/assets/112875495/94274d5d-1384-43bf-b981-10158d3235af">

`Priors`
<img width="1064" alt="image" src="https://github.com/LKFrederiksen/Cryptocarya/assets/112875495/3d9d3411-2ab2-4727-a082-5d04b027f1e7">

`Operators`
<img width="1067" alt="image" src="https://github.com/LKFrederiksen/Cryptocarya/assets/112875495/d9731290-16d7-49db-87f0-f69f90ed2111">

`MCMC`
<img width="1069" alt="image" src="https://github.com/LKFrederiksen/Cryptocarya/assets/112875495/bd395850-4f2b-4a91-9759-ec501caa5ed3">

See [`xml`]()

## Running BEAST
Molecular dating was performed using BEAST v1.10.4 (Suchard et al., 2018).

Tracer v1.7.2 was used to check for convergence of the runs (Rambaut et al., 2018) after discarding the initial 10% as burn-in. The analysis was run twice for up to 500 million generations until the ESS of all parameters was >200, logging every 10.000 generations after which log and tree files were combined in LogCombiner (Drummond et al., 2012). The median node ages with 95% highest posterior density (HPD) intervals were summarised and a maximum clade credibility tree was annotated with TreeAnnotator v1.10.4 (Drummond et al., 2012) and visualized using FigTree v1.4.4 (Rambaut, 2009).
