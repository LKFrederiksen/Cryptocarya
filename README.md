# **Cryptocarya group**

This repository contains all the scripts, including workflows I used to conduct my analyses, during my MSc thesis at Aarhus University. Additionally, the supplementary information for an article manuscript. The submitted MSc thesis can be read [here]().

During my thesis I studied the Cryptocarya group (Lauraceae), a pantropical group mainly linked to tropical rainforest (TRF), consisting of probably >3500 species. The focus of my thesis was understanding where the Cryptocarya group most likely originated and when, where and how this group achieved its pantropical distribution.
The thesis is structured in two parts: 
- Part A is a general introduction of the subjects treated in this thesis, which does not include any scripts.
- Part B is the article manuscripy entitled 'Long-distance dispersal shaped the pantropical distribution of the geographically structured Cryptocarya group', where I generated a phylogeny of the Cryptocarya group, dated the phylogeny using nine fossil calibrations and conducted ancestral range estimation.
    - The supplementary information can be found in the folder [`/Supplementary_information`](https://github.com/LKFrederiksen/Cryptocarya/tree/c75be30004c57bcb3f6e78a846e54deeddb694b0/Supplementary_information).
    - See below for a description of the workflows used to conduct the analyses included in the manuscript:

## **PART B: Phylogenomics, molecular dating and biogeography of the Cryptocarya group**

Workflows are found in the following folders:
     
- [1_Phylogeny](https://github.com/LKFrederiksen/Cryptocarya/tree/6f4d16c606cc2b1afc21b40ea2439306688dddc0/1_Phylogeny): Workflow to reproduce generated ASTRAL species tree of the Cryptocarya group, using the Angiosperms353 probe set [(Johnson et al., 2019)](https://www.ncbi.nlm.nih.gov/pubmed/30535394).
     
- [2_Molecular_dating](https://github.com/LKFrederiksen/Cryptocarya/tree/6f4d16c606cc2b1afc21b40ea2439306688dddc0/2_Molecular_dating)): Workflow to date the species tree using BEAST [(Drummond et al., 2018)](https://www.ncbi.nlm.nih.gov/pubmed/29942656).
      
- [3_Biogeography](https://github.com/LKFrederiksen/Cryptocarya/tree/6f4d16c606cc2b1afc21b40ea2439306688dddc0/3_Biogeography): Workflow to undertake ancestral range estimation and biogeograhical stochastical mapping (BSM) using BioGeoBEARS [(Matzke, 2013)](https://escholarship.org/content/qt44j7n141/qt44j7n141_noSplash_a25fcfcd2e4c86c599d6f7da5ee1d7f8.pdf?t=pga0we).

## Miscellaneous

The folder [`\Miscellaneous`](https://github.com/LKFrederiksen/Cryptocarya/tree/142f6f032675b593ef73c2ea7d9b6512ece5f551/Miscellaneous) contains various scripts that may be helpful. E.g., for visualisation of phylogenetic trees or to generate world maps at different points in time.

## For your information

Please do not expect all the scipts to run out of the box from your device.
