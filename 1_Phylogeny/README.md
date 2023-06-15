# Cryptocarya group phylogeny

## Getting ready - INSTALLATION 
In order to be able to run the pipeline, you first need to install the following programs:

I installed the packages via conda in separate environments.

### Programs 

FastQC version 0.11.9 (Andrews, 2010; )

MultiQC version  1.13 (Ewels et al., 2016; )

Hybpiper version 2.0.1 (Johnson et al. 2016; https://github.com/mossmatters/HybPiper)

BWA version 0.7.17 (Li & Durbin 2009; https://github.com/lh3/bwa)

Samtools version 1.7 (Li et al., 2009; https://github.com/samtools/samtools)

Picard version 2.18.29 (Picard Tool Kit, 2019; https://github.com/broadinstitute/picard)

MAFFT version 7.508 (Katoh & Standley, 2013; https://gitlab.com/sysimm/mafft)

Trimal version 1.4.1 (Capella-Gutiérrez et al. 2009; https://github.com/inab/trimal)

CIalign version 1.0.18 (Tumescheit et al. 2022; https://github.com/KatyBrown/CIAlign)

TAPER version 1.8.2 (Zhang et al. 2021; https://github.com/chaoszhang/TAPER)

IQtree version 2.2.0.3 (Nguyen et al. 2015; https://github.com/Cibiv/IQ-TREE)

ASTRAL-III version 3.2 (Yin et al. 2019; https://github.com/smirarab/ASTRAL/tree/MP)

Rstudio version 4.1.3 (RStudio Team 2023; http://www.rstudio.com/)

I ran the entire workflow via the workflow manager [GWF](https://github.com/gwforg/gwf)

### Directories
On my cluster account,on GenomeDK I ran all directories from the following path: /Workflow/Final_tree/

RAW data was put in the following two folders:
"/home/laurakf/cryptocarya/RawData/Lauraceae # Samples from Professor Jens G. Rohwer, Hamburg University, prepared for sequencing at Aarhus University.
- This data is available from [https://www.ncbi.nlm.nih.gov/sra/PRJNA939499](https://www.ncbi.nlm.nih.gov/sra/PRJNA939499), upon release.
"/home/laurakf/cryptocarya/RawData/PAFTOL # Data from [Baker et al., 2022](https://www.ncbi.nlm.nih.gov/pubmed/33983440).

### Target file
I used the Mega353 target file: [(https://github.com/chrisjackson-pellicle/NewTargets#newtargets)](https://github.com/chrisjackson-pellicle/NewTargets#newtargets)

## Running the pipeline

### The entire pipeline can be found in the following scripts [](). All paths and settings are available.

### Workflow overview

  1. Copy sequences from raw data folder to working directory. 
  2. FastQC quality check of the sequences with slidingwindow or (maxinfo)
  3. MultiQC wrapper to combine FastQC reports in one html file
  4. Trimming of the sequences using trimmomatic
     - Used following setting to trim the raw sequences.
    trimmomatic PE -threads 16 -phred33 {path_in}_R1.fastq {path_in}_R2.fastq -baseout {output}.fastq\
    ILLUMINACLIP:/home/laurakf/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:1:30:7:1:true\
    LEADING:30\
    SLIDINGWINDOW:4:30\
    MINLEN:40\
  6. FastQC quality check of the trimmed sequences
  7. MultiQC wrapper to combine FastQC reports in one html file
  8. Move files
     At this step all the 'PAFTOL' outgroup samples were moved to the HybPiper folder.
      - The remaining samples 'ingroup' were moved to the folder

  9. Hybpiper, using assemble (including intronerate) to get exons and introns (supercontigs) for each sample.
        Retrieving stats.
        Paralog retriever to check for paralogs.
      
  10. Coverage in 2 steps:
      1. Trimming for Coverage of sequencing and joining the exon and intron of each species to create supercontigs.
         Here I used the script `coverage_Laura-py`
      2. Trimming for Coverage of sequencing and joining the exon and intron to create supercontigs for the 10% best genes and creating exons for the 90% 'worst' genes as determined by an earlier conducted sortadate run.
          Here I used the script `coverage_comp.py`
 
  11. Using script sample2genes (pastataster) to retrieve the gene files including all sample sequences found for each gene.
  
  12. MAFFT to align sequences
      
  14.  Trimal to trim alignments

  15. AMAS to get summary statistics to trim alignments.
  
  16. Optrimal to trim alignments
  
  17. Move Optrimal files
 
  18. CIalign to trim alignments
  
  19. Taper to trim alignments
  
  20. Move HybPiper outgroup data into same folder as ingroup
  
  21. Exon mapper # Map best exons to alignment to enable using PartitionFinder in IQtree
  
  22. Copying alignments and makin gpartitions
  
  23. IQtree to generate a hypothesis for each gene and generate gene trees
  
  24. Generate ASTRAL tree from gene trees

## Contig selection

I furthermore followed the following protocol to deter mine, which genes were best.
In our pipeline to generate a species tree of the Cryptocarya group including outgroup we first generated a tree, using supercontigs for everything (ingroup and outgroup). However, we got an unrealistic topology in the outgroup contradicting Baker et al. (2022) and additional studies (Chanderbali et al., 2001; Song et al., 2019), probably attributable to alignment issues in the non-exonic regions. We used Baker et al. (2022) as benchmark to remove non-exonic regions. Rerunning the analysis with exons we found that the non-exonic regions caused this issue. As compromise we generated a species tree from supercontigs but removing non-exonic regions from the outgroup as necessary. We did this by using SortaDate (Smith et al., 2018) to first identify which of the supercontig gene trees that were most like the outgroup topology of a previously generated exon species tree (only exons in ingroup and outgroup). Hereafter, we generated new gene trees (supercontig + exon) for the outgroup by running IQ-TREE (Minh et al., 2020) and combining different combinations of the previously generated supercontigs and exons until the anticipated outgroup topology was achieved. The correct outgroup included the 10% genes with the highest bipartition scores as supercontig genes and the remaining 90% as exons (Appendix S2). Finally, we reran the analyses using supercontigs for the ingroup and the abovementioned combination of supercontigs and exons for the outgroup.  

