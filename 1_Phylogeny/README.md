# Cryptocarya group phylogeny

# Getting ready - INSTALLATION 


#----------------------------------------------------------------------------------------------------------------------------
# ### Workflow 1 ###
# ----------------------------------------------------------------------------------------------------------------------------
# 1: FastQC quality check of the sequences with slidingwindow or (maxinfo)
# 2: MultiQC wrapper to combine FastQC reports in one html file
# 3: Trimming of the sequences using trimmomatic
# 4: FastQC quality check of the trimmed sequences
# 5: MultiQC wrapper to combine FastQC reports in one html file
# 6: Hybpiper, using assemble (including intronerate) to get exons and introns for each sample. Retrieving stats. Paralog retriever to check for paralogs.
# 7: Coverage in 2 steps: 
# 	- 1. Trimming for Coverage of sequencing and joining the exon and intron of each species to create supercontigs.
#   - 2. Trimming for Coverage of sequencing and joining the exon and intron to create supercontigs for the 10% best genes and creating exons for the 90% 'worst' genes as determined by an earlier conducted sortadate run.
# 8: Using script sample2genes (pastataster) to retrieve the gene files including all sample sequences found for each gene.
# 9: MAFFT to align sequences
# 10: Trimal to trim alignments
# 11: Optrimal to trim alignments
# 12: CIalign to trim alignments
# 13: Taper to trim alignments
# 14: IQtree to generate a hypothesis for each gene and generate gene trees
# 15: Generate ASTRAL tree from gene trees


Filer skal køres fra et directory der hedder følgende: /Workflow/Final_tree/

"/home/laurakf/cryptocarya/RawData/PAFTOL/"
