# Dating the phylogeny with BEAST

[BEAST](https://beast.community) [(Drummond & Rambaut, 2007)](https://www.ncbi.nlm.nih.gov/pubmed/17996036
https://bmcecolevol.biomedcentral.com/counter/pdf/10.1186/1471-2148-7-214.pdf) based on Bayesian inference was used to date the phylogenetic tree. To do this we used the three best nuclear genes, nine fossils (primary calibrations) and one secondary calibration on the root (see Table ... in [MSc thesis]()). Moreover, we fixed the tree topology to that of the ASTRAL tree.

Look at part 3 of my [`workflow.py`]() file based on [gwf]() to see the steps and code I used to undertake parts of the molecular dating analysis in the [GenomeDK](https://genome.au.dk/docs/) cluster. 

## Selecting fossils
We selected the fossils based on [Ramirez-Barahona et al., 2022](https://www.nature.com/articles/s41559-020-1241-3) and the fossils included in [Li et al. (2020)](https://www.ncbi.nlm.nih.gov/pubmed/32619568), except the fossil [_Cryptocaryoxylon gippslandicum_](https://doi.org/10.1080/03115518608619157). See [`/Fossil_table.pdf`](https://github.com/LKFrederiksen/Cryptocarya/blob/ffce9223b5dadbb4187c49fbbc932fc810e2b024/2_Molecular_dating/Fossil_table.pdf).

## Selecting genes
To select genes for dating I used [SortaData](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0197433&type=printable) to get a list of the 'best' genes, thus sorted after the genes with least root to tip variance and best bipartition support compared to the species tree. See [`workflow.py`](). This is done in order to select the genes, whose gene trees have the least topological conflict with the species tree. Here the three genes (5594, 5620 and 6139) were selected.

## Alignments
Based on the output of SortaDate I selected the three 'best' genes that contained 'supercontigs' (exons including the flanking intron and spacer-regions) in both the ingroup and outgroup. 
For each gene in the alignments, if a sample dit not have a sequence, as no reads had been recovered for the specific gene, then gaps were added as required for BEAST to run. 

This was done using [`/add_seqs.zip`](2_Molecular_dating/add_seqs.zip) based on Biopython v1.7.5, an adapted version of a script [/1.1_add_empty_seqs_for_sp.py](2_phylo_dating/1.1_add_empty_seqs_for_sp.py)created by Maya Schrödl.
This folder should be unzipped and the python script `add_empty_seqs_for_sp.py` should be run from within the folder `add_seqs`.

For the species where more than one sample was available, only one individual was chosen to comply with BioGeoBEARS. The individual with the alignment that had the least gaps was selected.
These genes were manually cleaned using the program [AliView](https://doi.org/10.1093/bioinformatics/btu531). 

The cleaned alignments can be retrieved from the folder `alignments_cleanTrim` in the [`BEAUti.zip`](2_Molecular_dating/BEAUti.zip) folder.

The tips of the samples that were in excess were dropped from the ASTRAL species tree, so it could be used as starting tree and to fix the topology.

## Preparing the xml file in BEAUti
All files used in BEAUti can be found in the folder [`2_Molecular_dating/BEAUti.zip`](2_Molecular_dating/BEAUti.zip).

To run the molecular dating analysis in BEAST a xml file is required. This file can be generated in the BEAST GUi [BEAUti](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3408070/pdf/mss075.pdf). I used version 1.10.4. 

The best-performing evolutionary model for each gene was determined using the Akaike Information Criterion (AIC; Akaike, 1974) as implemented in IQ-TREE ModelFinder [(Kalyaanamoorthy et al., 2017)](https://www.nature.com/articles/nmeth.4285.pdf). The GTR+Γ+I with 4 Γ categories was selected for gene 5620 and 5594, whereas HKY+ Γ with 4 Γ categories was selected for gene 6139. 

An uncorrelated relaxed clock model assuming an uncorrelated lognormal rate prior and a birth-death speciation prior was selected. The [`pruned ASTRAL tree`](2_Molecular_dating/BEAUti.zip) was selected as user-specified starting tree  and used as constraint by deleting the commands: subtreeSlide; narrowExchange; wideExchange; wilsonBalding in BEAUti v.1.10.4, because a topology inferred from 300+ genes is more reliable than one inferred from three.

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

See `3_bestGenes` in [`/BEAUti.zip`](https://github.com/LKFrederiksen/Cryptocarya/blob/ffce9223b5dadbb4187c49fbbc932fc810e2b024/2_Molecular_dating/BEAUti.zip).

## Running BEAST
Input files for BEAST can be found in the folder [`/BEAST_input.zip`](https://github.com/LKFrederiksen/Cryptocarya/blob/ffce9223b5dadbb4187c49fbbc932fc810e2b024/2_Molecular_dating/BEAST_input.zip).

Molecular dating was performed using BEAST v1.10.4 [(Suchard et al., 2018)](https://watermark.silverchair.com/vey016.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAArkwggK1BgkqhkiG9w0BBwagggKmMIICogIBADCCApsGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMDF0SPJ1oEqdfjGEfAgEQgIICbMqhzaqUrxMpYR_pw8Ql0OX8-Y52fxooxSOSIUprGj3LXc5gMWwkdE-yM6rxDQM7nPMElxYJ0I7w8dxgxIHTlFzUkblXAlM6qzvrlfeL5Yr1Fqv9m74Pm7kBU1kQCd1lJy9vt3GoWe16uAuCSdXfNumW4O4QT3pPDR-xu8RcZ-eEoQuQxLORjOlP4QoXelfoc0gWPLSI-O7dXYBEFHV7_OTGbNN5Mx2i7BcEuElSu0iOBl8iobeFQ60Vmj-E-sYpw3McFNLHzpVcbuiH0DXuOZX9tcqNUS5iPWyZPn2-KTZEhi6OxEGxiPrZzeYs9OvmFzSm42LhqEp7c_Nuf3VdDWZIoxlzltmkWzU7s62AK2KclnnvgnEBT8B5uEVnwkgFx9Q3Iy8mwjzkpV35mx9s0KssTO2cEkIrXIHyLv6_7vot10bL5SzhUAfU4NG2imx_tsrnGtdqTOLxs4MC-1543oBhgCCY8XGWZyGaMwjPQfklVKL1BgYLUi9mWWHKoxLNfEArjt7xnBNyV0ck2AT1iQmnfd1wkMOEkveFZa0QaXNB0HDjnFt_xr7Rh23fUXo-C--EH8HfGVKhP8u2OR1Z02hlTS7ZYLLfr5GYdBeVUqZNtfe-3N3m3HZWFFFiK5Ud0KcA5htQ6TG-3l2HfRrfmsJkBaYLtGZra9raj5aVzg_hlbxmumFbGYaW9HkzXmgzZxTvd2k7_7Qsf5FmKIueu_bpil2sy1B1-Wc-3C-bKQbVqWeBDDWGUZoLRviDhDrqqfkHQKvQ_Y7sNgbeOL7_fNo-mbBRmVQnjTYGijdn0eWVNLFJBHqIDhaNepyk). I ran BEAST through the GenomeDK cluster. See part 3 of [`workflow.py`]() for an example of how the program can be run using the xml files.

Tracer v1.7.2 [(Rambaut et al., 2018)](https://www.ncbi.nlm.nih.gov/pubmed/29718447) was used to check for convergence of the runs  after discarding the initial 10% as burn-in. The analysis was run twice for up to 500 million generations until the ESS of all parameters was >200, logging every 10.000 generations after which log and tree files were combined in LogCombiner [(Drummond et al., 2012)](). The median node ages with 95% highest posterior density (HPD) intervals were summarised and a maximum clade credibility tree was annotated with TreeAnnotator v1.10.4 [(Drummond et al., 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22367748
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3408070/pdf/mss075.pdf) and visualized using FigTree v1.4.4 [(Rambaut, 2009)](http://tree.bio.ed.ac.uk/software/figtree/).
