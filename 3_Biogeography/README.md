# **Biogeographical analyses with BioGeoBEARS**

All biogeographical analyses were undertaken in the programming language [R]() using  [BioGeoBears](https://github.com/nmatzke/BioGeoBEARS) v. [(Matzke et al., 2013)](). Documentation for all the methods used can be found [here](http://phylo.wikidot.com/biogeobears). The outgroup was excluded from these analyses.

For each of the three analyses a folder including all required files to run the analyses has been uploaded. To run the analyses, the analyses should be run within each of the respective folders. Remember to change paths in the script within the folder to the location of the folder on your device. 

For the 'Ancestral range estimation' the folders also include the output. This is not the case for the BSM analysis, as the output is too space consuming to be uploaded.

## Ancestral range estimation

We defined phytogeographical regions based on [Carta et al. (2022)](https://doi.org/10.1111%2Fnph.17844) and assigned extant species based on the [World Checklist of Vascular Plants](https://checklistbuilder.science.kew.org/reportbuilder.do) [(Govaerts et al., 2021)](https://www.nature.com/articles/s41597-021-00997-6.pdf).

Three maximum likelihood (ML) models of ancestral range evolution (DEC, DIVALIKE; BAYAREALIKE) were applied in this study, also including the 'jump dispersal' parameter _j_. Additionally, we conducted time-stratified (TS) analyses by splitting the tree into three time slices (TS1: 80-40 Ma, TS2: 40-20 Ma, TS1: 20-0 Ma) to assign different dispersal rate scalers to each time slice. This was done to account for varying dispersal probabilities through time. The dispersal rate scalers were assigned based on the rules set up in [(Touissant et al., 2017)](https://doi.org/10.1111/jbi.12977). In total 12 models were tested and compared using AIC [(Burnham & Anderson, 2002)](https://doi.org/10.1007/b97636) to select the most likely model.

### [`3_Biogeography/10_GeoAreas_3genes.zip`](3_Biogeography/10_GeoAreas_3genes.zip)

The script `BioGeoBEARS_3genes.R` within the folder is based on the [example BioGeoBEARS script](http://phylo.wikidot.com/biogeobears#script) and runs all the analyses without time stratification and dispersal rate scalers.

Thus, `BioGeoBEARS_3genes.R` within the folder runs the following six models:
- DEC
- DEC + _j_
- DIVALIKE
- DIVALIKE + _j_
- BAYAREALike
- BAYAREALIKE + _j_`

### [`3_Biogeography/10_GeoAreas_TS_3genes.zip`](3_Biogeography/10_GeoAreas_TS_3genes.zip)

This script is the same as `BioGeoBEARS_TS_3genes.R`, but is time-stratified and includes dispersal rate scalers.

Thus, `BioGeoBEARS_TS_3genes.R` within this folder runs the following six models:
- DEC + TS
- DEC + TS + _j_
- DIVALIKE + TS
- DIVALIKE + TS + _j_
- BAYAREALike + TS
- BAYAREALIKE + TS + _j_`

## Biogeographical Stochastic Mapping (BSM)

BSM is programmed to do stochastic mapping on any of the abovementioned biogeographical models (and more variables). In this study BSM was used to explore which speciation mechanisms were dominant in the different time slices and the number of dispersal events (range expansion, d, or jump dispersal, _j_) happening within the different time slices and between the defined biogeographical regions.

### [`BSM individual slice script`]

The BSM was therefore conducted for individual time slices using the script `BioGeoBEARS_TS_BSM_3genes.R` based on the following [example script](http://phylo.wikidot.com/biogeographical-stochastic-mapping-example-script#BSM_script).
