
# ----------------------------------------------------------------------------------------------------------------------------
# This workflow is used to transform the raw sequence data into sequences ready for alignment.
# Workflow is the following

# 1: Secapr quality check of the sequences
# 2: Trimming of the sequences using trimmomatic
# 3: Hybpipering the species in order to create the exons of the baits for each species
# 4: Checking for Paralogs
# 5: Running Intronerate again to get the Introns for each species
# 6: Trimming for Coverage of sequencing and joining the exon and intron of each species to create supercontigs

# ----------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------------
# This code is a variant of species_workflow.py by Oscar Balslev Wrisberg and workflow_completed.py by Paola de Lima Ferreira 
# Edited by Laura Kragh Frederiksen 12/10/2022
# ----------------------------------------------------------------------------------------------------------------------------

from os import O_SYNC, name
from gwf import Workflow
import os.path
# import math
# import glob

gwf = Workflow()

########################################################################################################################
###############################################---- Fastqc quality check raw ----#######################################
########################################################################################################################
def fastqc_raw(name,path_in ,path_out, done,):
    """Quality checking using fastqc as this should work on individual species"""
    path_ins = []
    outputs = [path_out+name+"_R1_fastqc.html",path_out+name+"_R2_fastqc.html", done]
    options = {'cores': 1, 'memory': "8g", 'walltime': "00:30:00", 'account':"cryptocarya"}


    spec = """

    echo {name}

    source /home/laurakf/miniconda3/etc/profile.d/conda.sh

    conda activate fastqc

    fastqc -o {path_out} {path_in}{name}_R1.fastq {path_in}{name}_R2.fastq
    
    echo touching {done}
    touch {done}

    """.format(path_in = path_in, name = name, path_out = path_out, done = done)

    return (path_ins, outputs, options, spec)

########################################################################################################################
##############################################---- Multiqc quality check raw ----#######################################
########################################################################################################################
def multiqc_raw(path_in ,path_out, done,):
    """Quality checking using multiqc"""
    inputs = [path_in+"Ocotea-foetens-WE521_R1_fastqc.html",path_in+"Ocotea-gabonensis-WE522_R1_fastqc.html",path_in+"Ocotea-meziana-WE523_R1_fastqc.html",path_in+"Pleurothyrium-cuneifolium-WE524_R1_fastqc.html",path_in+"Mespilodaphne-cymbarum-WE525_R1_fastqc.html",path_in+"Damburneya-gentlei-WE526_R1_fastqc.html",path_in+"Ocotea-glaucosericea-WE527_R1_fastqc.html",path_in+"Ocotea-complicata-WE528_R1_fastqc.html",path_in+"Ocotea-javitensis-WE529_R1_fastqc.html",path_in+"Ocotea-skutchii-WE530_R1_fastqc.html",path_in+"Ocotea-sinuata-WE531_R1_fastqc.html",path_in+"Ocotea-botrantha-WE532_R1_fastqc.html",path_in+"Nectandra-lineatifolia-WE533_R1_fastqc.html",path_in+"Ocotea-foetens-WE521_R2_fastqc.html",path_in+"Ocotea-gabonensis-WE522_R2_fastqc.html",path_in+"Ocotea-meziana-WE523_R2_fastqc.html",path_in+"Pleurothyrium-cuneifolium-WE524_R2_fastqc.html",path_in+"Mespilodaphne-cymbarum-WE525_R2_fastqc.html",path_in+"Damburneya-gentlei-WE526_R2_fastqc.html",path_in+"Ocotea-glaucosericea-WE527_R2_fastqc.html",path_in+"Ocotea-complicata-WE528_R2_fastqc.html",path_in+"Ocotea-javitensis-WE529_R2_fastqc.html",path_in+"Ocotea-skutchii-WE530_R2_fastqc.html",path_in+"Ocotea-sinuata-WE531_R2_fastqc.html",path_in+"Ocotea-botrantha-WE532_R2_fastqc.html",path_in+"Nectandra-lineatifolia-WE533_R2_fastqc.html"] 
    outputs = [path_out+"multiqc_report.html", done]
    options = {'cores': 1, 'memory': "8g", 'walltime': "00:30:00", 'account':"cryptocarya"}


    spec = """

    source /home/laurakf/miniconda3/etc/profile.d/conda.sh

    conda activate multiqc

    multiqc -o {path_out} {path_in}
    
    touch {done}

    """.format(path_in = path_in, path_out = path_out, done = done)

    return (inputs, outputs, options, spec)


#######################################################################################################################################
################################################---- Trimmomatic - SLidingwindow----###################################################
#######################################################################################################################################
def trimmomatic(name, path_in, path_out, done):
    """Trimming raw data using trimmomatic with custom adapter.
    Afterwards combines paired and unpaired reads for forward and reverse reads respectively for each species 
    to enable post-trimming secapr quality_check for comparability before and after trimming """
    path_ins = []
    outputs = [path_out+name+"_UN.fastq",path_out+"secapr_postrim/"+name+"_UN.fastq", done, path_out+name+"_1P.fastq", path_out+name+"_2P.fastq"]
    options = {'cores': 16, 'memory': "10g", 'walltime': "01:00:00", 'account':"cryptocarya"}

    spec = """
    source /home/laurakf/miniconda3/etc/profile.d/conda.sh
    conda activate trimmomatic

    trimmomatic PE -threads 16 -phred33 {path_in}_R1.fastq {path_in}_R2.fastq -baseout {output}.fastq\
    ILLUMINACLIP:/home/laurakf/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:1:30:7:1:true\
    LEADING:30\
    SLIDINGWINDOW:4:30\
    MINLEN:40\
    2>> stderr_trim_loop_output.txt

    echo combining {path_out}{name}_1P.fastq and {path_out}{name}_1U.fastq into {path_out}secapr_postrim/{name}_1PU.fastq 
    cat {path_out}{name}_1P.fastq {path_out}{name}_1U.fastq > {path_out}secapr_postrim/{name}_1PU.fastq 

    echo combining {path_out}{name}_2P.fastq and {path_out}{name}_2U.fastq into {path_out}secapr_postrim/{name}_2PU.fastq 
    cat {path_out}{name}_2P.fastq {path_out}{name}_2U.fastq > {path_out}secapr_postrim/{name}_2PU.fastq

    echo combining {path_out}{name}_1U.fastq {path_out}{name}_2U.fastq > {path_out}{name}_UN.fastq
    cat {path_out}{name}_1U.fastq {path_out}{name}_2U.fastq > {path_out}{name}_UN.fastq
    cp {path_out}{name}_UN.fastq {path_out}secapr_postrim/


    echo Removing {path_out}{name}_1U.fastq
    rm {path_out}{name}_1U.fastq

    echo Removing {path_out}{name}_2U.fastq
    rm {path_out}{name}_2U.fastq

    touch {done}

    """.format(path_in = path_in+name, output = path_out+name, done = done, name = name, path_out = path_out)

    return (path_ins, outputs, options, spec)


#######################################################################################################################################
################################################---- Trimmomatic - Maxinfo ----########################################################
#######################################################################################################################################
# def trimmomatic(name, path_in, path_out, done):
#     """Trimming raw data using trimmomatic with custom adapter.
#     Afterwards combines paired and unpaired reads for forward and reverse reads respectively for each species 
#     to enable post-trimming secapr quality_check for comparability before and after trimming """
#     path_ins = []
#     outputs = [path_out+name+"_UN.fastq",path_out+"secapr_postrim/"+name+"_UN.fastq", done, path_out+name+"_1P.fastq", path_out+name+"_2P.fastq"]
#     options = {'cores': 16, 'memory': "10g", 'walltime': "01:00:00", 'account':"cryptocarya"}

#     spec = """
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate trimmomatic

#     trimmomatic PE -threads 16 -phred33 {path_in}_R1.fastq {path_in}_R2.fastq -baseout {output}.fastq\
#     ILLUMINACLIP:/home/laurakf/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:1:30:7:1:true\
#     LEADING:3\
#     TRAILING:3\
#     MAXINFO:40:0.8\
#     MINLEN:36\
#     2>> stderr_trim_loop_output.txt

#     echo combining {path_out}{name}_1P.fastq and {path_out}{name}_1U.fastq into {path_out}secapr_postrim/{name}_1PU.fastq 
#     cat {path_out}{name}_1P.fastq {path_out}{name}_1U.fastq > {path_out}secapr_postrim/{name}_1PU.fastq 

#     echo combining {path_out}{name}_2P.fastq and {path_out}{name}_2U.fastq into {path_out}secapr_postrim/{name}_2PU.fastq 
#     cat {path_out}{name}_2P.fastq {path_out}{name}_2U.fastq > {path_out}secapr_postrim/{name}_2PU.fastq

#     echo combining {path_out}{name}_1U.fastq {path_out}{name}_2U.fastq > {path_out}{name}_UN.fastq
#     cat {path_out}{name}_1U.fastq {path_out}{name}_2U.fastq > {path_out}{name}_UN.fastq
#     cp {path_out}{name}_UN.fastq {path_out}secapr_postrim/


#     echo Removing {path_out}{name}_1U.fastq
#     rm {path_out}{name}_1U.fastq

#     echo Removing {path_out}{name}_2U.fastq
#     rm {path_out}{name}_2U.fastq

#     touch {done}

#     """.format(path_in = path_in+name, output = path_out+name, done = done, name = name, path_out = path_out)

#     return (path_ins, outputs, options, spec)


########################################################################################################################
###########################################---- Fastqc quality check trimmed ----#######################################
########################################################################################################################
def fastqc_trimmed(name,path_in ,path_out, done,):
     """Quality checking using fastqc as this should work on individual species"""
     path_ins = [path_in+name+"_UN.fastq", path_in+name+"_1PU.fastq", path_in+name+"_2PU.fastq"] # The files gwf looks for before it runs.
     outputs = [path_out+name+"_1PU_fastqc.html", path_out+name+"_2PU_fastqc.html",path_out+name+"_UN_fastqc.html", done]
     options = {'cores': 1, 'memory': "8g", 'walltime': "00:30:00", 'account':"cryptocarya"}


     spec = """

     echo {name}
     
     source /home/laurakf/miniconda3/etc/profile.d/conda.sh

     conda activate fastqc

     fastqc -o {output} {path_in}{name}_1PU.fastq {path_in}{name}_2PU.fastq {path_in}{name}_UN.fastq
    
     touch {done}

     """.format(path_in = path_in, name = name, output = path_out, done = done)

     return (path_ins, outputs, options, spec)


# ########################################################################################################################
# ##########################################---- Multiqc quality check trimmed ----#######################################
# ########################################################################################################################
# def multiqc_trimmed(path_in ,path_out, done,):
#     """Quality checking using multiqc"""
#     inputs = [path_in+"Ocotea-foetens-WE521_1PU_fastqc.html",path_in+"Ocotea-gabonensis-WE522_1PU_fastqc.html",path_in+"Ocotea-meziana-WE523_1PU_fastqc.html",path_in+"Pleurothyrium-cuneifolium-WE524_1PU_fastqc.html",path_in+"Mespilodaphne-cymbarum-WE525_1PU_fastqc.html",path_in+"Damburneya-gentlei-WE526_1PU_fastqc.html",path_in+"Ocotea-glaucosericea-WE527_1PU_fastqc.html",path_in+"Ocotea-complicata-WE528_1PU_fastqc.html",path_in+"Ocotea-javitensis-WE529_1PU_fastqc.html",path_in+"Ocotea-skutchii-WE530_1PU_fastqc.html",path_in+"Ocotea-sinuata-WE531_1PU_fastqc.html",path_in+"Ocotea-botrantha-WE532_1PU_fastqc.html",path_in+"Nectandra-lineatifolia-WE533_1PU_fastqc.html",path_in+"Ocotea-foetens-WE521_2PU_fastqc.html",path_in+"Ocotea-gabonensis-WE522_2PU_fastqc.html",path_in+"Ocotea-meziana-WE523_2PU_fastqc.html",path_in+"Pleurothyrium-cuneifolium-WE524_2PU_fastqc.html",path_in+"Mespilodaphne-cymbarum-WE525_2PU_fastqc.html",path_in+"Damburneya-gentlei-WE526_2PU_fastqc.html",path_in+"Ocotea-glaucosericea-WE527_2PU_fastqc.html",path_in+"Ocotea-complicata-WE528_2PU_fastqc.html",path_in+"Ocotea-javitensis-WE529_2PU_fastqc.html",path_in+"Ocotea-skutchii-WE530_2PU_fastqc.html",path_in+"Ocotea-sinuata-WE531_2PU_fastqc.html",path_in+"Ocotea-botrantha-WE532_2PU_fastqc.html",path_in+"Nectandra-lineatifolia-WE533_2PU_fastqc.html",path_in+"Ocotea-foetens-WE521_UN_fastqc.html",path_in+"Ocotea-gabonensis-WE522_UN_fastqc.html",path_in+"Ocotea-meziana-WE523_UN_fastqc.html",path_in+"Pleurothyrium-cuneifolium-WE524_UN_fastqc.html",path_in+"Mespilodaphne-cymbarum-WE525_UN_fastqc.html",path_in+"Damburneya-gentlei-WE526_UN_fastqc.html",path_in+"Ocotea-glaucosericea-WE527_UN_fastqc.html",path_in+"Ocotea-complicata-WE528_UN_fastqc.html",path_in+"Ocotea-javitensis-WE529_UN_fastqc.html",path_in+"Ocotea-skutchii-WE530_UN_fastqc.html",path_in+"Ocotea-sinuata-WE531_UN_fastqc.html",path_in+"Ocotea-botrantha-WE532_UN_fastqc.html",path_in+"Nectandra-lineatifolia-WE533_UN_fastqc.html"] 
#     outputs = [path_out+"multiqc_report.html", done]
#     options = {'cores': 1, 'memory': "8g", 'walltime': "00:30:00", 'account':"cryptocarya"}


#     spec = """

#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh

#     conda activate multiqc

#     multiqc -o {path_out} {path_in}
    
#     touch {done}

#     """.format(path_in = path_in, path_out = path_out, done = done)

#     return (inputs, outputs, options, spec)



########################################################################################################################
################################################---- Hybpiper ----######################################################
########################################################################################################################
def hybpiper(name, p1, p2, un, path_out, path_in, done):
    """Hybpiper."""
    path_ins = [path_in+name+p1, path_in+name+p2, path_in+name+un] # The files which the job will look for before it runs
    outputs = [path_out+name, done] # The files which will have to be created in order for the job to be "completed"
    options = {'cores': 2, 'memory': "8g", 'walltime': "10:00:00", 'account':"cryptocarya"} #Slurm commands

    spec = """

    source /home/laurakf/miniconda3/etc/profile.d/conda.sh

    conda activate HybPiper

    TMPDIR=/scratch/$SLURM_JOBID
    export TMPDIR
    mkdir -p $TMPDIR
    cd $TMPDIR
    
    hybpiper assemble --cpu 2 --targetfile_dna /home/laurakf/cryptocarya/TargetFile/mega353.fasta --readfiles {p1} {p2} --unpaired {un} --prefix {name} --bwa --run_intronerate

    mv {name} /home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/

    touch {done}
    touch {out}{name}

    """.format(name=name, p1=path_in+name+p1, p2=path_in+name+p2, un=path_in+name+un, out=path_out+name, done=done)


    return (path_ins, outputs, options, spec)

########################################################################################################################
###################################################---- Stats ----######################################################
########################################################################################################################

# In this step you should run the statistics on the folder where we have the Hybpiper_results
# I did not create a folder just for Hybpiper results, then I will create here and move the assemble results to there

def stats(path_in, done, path_out, in_done, name):
   """Gather statistics about the HybPiper run(s).""" 
   path_ins = [path_in+name, in_done] # The files that has to be present before the job runs.
   outputs = [path_out+"seq_lengths.tsv", path_out+"hybpiper_stats.tsv", path_out+"recovery_heatmap.png"]  # The files which will have to be created in order for the job to be "completed"
   options = {'cores': 2, 'memory': "4g", 'walltime': "1:00:00", 'account':"cryptocarya"} #Slurm commands

   spec = """
   
   source /home/laurakf/miniconda3/etc/profile.d/conda.sh

   conda activate HybPiper
    
   cd {path_in}
    
   hybpiper stats --targetfile_dna /home/laurakf/cryptocarya/TargetFile/mega353.fasta 'gene' {path_in}namelist.txt # Get stats

   hybpiper recovery_heatmap {path_in}seq_lengths.tsv # Make heatmap

   mv seq_lengths.tsv {path_out} # Move all stats and the heatmap to a new subfolder
    
   mv hybpiper_stats.tsv {path_out}

   mv recovery_heatmap.png {path_out} 

   touch {done}
      
   """.format(path_in = path_in, done = done, path_out = path_out, in_done = in_done, name = name)

   return (path_ins, outputs, options, spec) 


# ########################################################################################################################
# #############################################---- Paralogs ----#########################################################
# ########################################################################################################################

def paralogs(name, path_in, done, in_done, path_out):
    """Run HybPiper v. 2.1 - paralog retriever """
    path_ins = [path_in+name, in_done]
    outputs = [done]
    options = {'cores': 2, 'memory': "10g", 'walltime': "0:30:00", 'account':"cryptocarya"}

    spec = """
    
    source /home/laurakf/miniconda3/etc/profile.d/conda.sh
    
    conda activate HybPiper
    
    cd {path_in}

    hybpiper paralog_retriever namelist.txt -t_dna /home/laurakf/cryptocarya/TargetFile/mega353.fasta
    
    mv paralog_report.tsv {path_out}
    mv paralogs_above_threshold_report.txt {path_out}
    mv paralogs_all {path_out}
    mv paralogs_no_chimeras {path_out}
    mv paralog_heatmap.png {path_out}

    touch {done}

     """.format(name = name, done = done, path_in = path_in, path_out = path_out, in_done = in_done)
    
    return (path_ins, outputs, options, spec)

# def paralogs(name, path_in, path_out, done, in_done):
#    """Find Paralog genes and write them on the file called paralog.txt"""
#    path_ins = [path_in+name, in_done]
#    outputs = [path_out+"paralog.txt", done]  # The files which will have to be created in order for the job to be "completed"
#    options = {'cores': 2, 'memory': "10g", 'walltime': "1:00:00", 'account':"cryptocarya"}

#    spec = """
   
#    source /home/laurakf/miniconda3/etc/profile.d/conda.sh

#    conda activate HybPiper

#    cd {path_in}

#    python /home/laurakf/cryptocarya/Programs/HybPiper-master/paralog_investigator.py {name} 2>> paralog.txt
   
#    mv paralog.txt {path_out}

#    touch {done}

#    """.format(name = name, done = done, path_in = path_in, path_out = path_out, in_done = in_done)

#    return (path_ins, outputs, options, spec)

def no_paralogs(name, path_in, done, no_paralogs, in_done):
    """Wrapper script to continue pipeline when Hybpiper finds no paralogs"""
    inputs = [path_in+name, in_done]
    outputs = [done]
    options = {'cores': 2, 'memory': "10g", 'walltime': "0:05:00", 'account':"cryptocarya"}

    spec = """

    touch {done}
    touch {np}

    """.format(done=done, np=no_paralogs, in_done = in_done)
    
    return(inputs, outputs, options, spec)



# ########################################################################################################################
# #############################################---- Coverage ----#########################################################
# ########################################################################################################################

#This script does the following:
# Gather all contigs from each sample in one fasta file: coverage/sample.fasta
# Map paired and unpaired reads to that fasta using BWA mem
# Deduplicate reads using Picard
# Calculate depth using samtools
# Mask/strip any bases with coverage <2
# Generate a new trimmed sample-level fasta: coverage/sample_trimmed.fasta

def coverage(name, path_in, path_out, done,all_bam,all_sorted_bam, all_sorted_bam_bai, bam, cov,fasta,fasta_amb,fasta_ann,fasta_bwt,fasta_pac,fasta_sa,trimmed_fasta,up_bam,dir_in,dir_out, dir_wrk):
    """Calculating coverage of sequences."""
    path_ins = [path_in+name]
    outputs = [path_out+name+all_bam,
     path_out+name+all_sorted_bam,
      path_out+name+all_sorted_bam_bai,
       path_out+name+bam,
    path_out+name+cov,
     path_out+name+fasta,
      path_out+name+fasta_amb,
       path_out+name+fasta_ann,
        path_out+name+fasta_bwt,
    path_out+name+fasta_pac,
     path_out+name+fasta_sa,
      path_out+name+trimmed_fasta,
       path_out+name+up_bam,done] #ALL the output files
    options = {'cores': 4, 'memory': "20g", 'walltime': "08:00:00", 'account':"cryptocarya"}

    spec = """
    
    source /home/laurakf/miniconda3/etc/profile.d/conda.sh
    
    conda activate HybPiper
    
    cd {path_in}

    python3 /home/laurakf/cryptocarya/Scripts/coverage.py {name} {dir_in} {dir_out} {dir_wrk}
    
    touch {done}

    """.format(name = name, done = done, path_in = path_in, dir_in = dir_in, dir_out = dir_out, dir_wrk = dir_wrk)

    return (path_ins, outputs, options, spec)

########################################################################################################################
#############################################---- Retrieve ----#########################################################
########################################################################################################################

#Think about doing blacklisting here? you could just remove species from the inputs here if you dont want them in the downstream analysis

def retrieve(path_in, done):
    """Retrieve gene sequences from all the species and create an unaligned multifasta for each gene."""
    path_ins = [path_in+"Pleurothyrium-cuneifolium-WE524_trimmed.fasta", path_in+"Ocotea-skutchii-WE530_trimmed.fasta", path_in+"Ocotea-sinuata-WE531_trimmed.fasta", path_in+"Ocotea-meziana-WE523_trimmed.fasta", path_in+"Ocotea-javitensis-WE529_trimmed.fasta", path_in+"Ocotea-glaucosericea-WE527_trimmed.fasta", path_in+"Ocotea-gabonensis-WE522_trimmed.fasta", path_in+"Ocotea-foetens-WE521_trimmed.fasta", path_in+"Ocotea-complicata-WE528_trimmed.fasta", path_in+"Mespilodaphne-cymbarum-WE525_trimmed.fasta", path_in+"Damburneya-gentlei-WE526_trimmed.fasta"]
    outputs = [done]
    options = {'cores': 4, 'memory': "5g", 'walltime': "1:00:00", 'account':"cryptocarya"}

    spec = """
    
    source /home/laurakf/miniconda3/etc/profile.d/conda.sh

    conda activate HybPiper

    cd {path_in}

    ls *trimmed.fasta > filelist.txt

    python3 /home/laurakf/cryptocarya/Scripts/sample2genes.py > outstats.csv

    touch {done}

    """.format(path_in = path_in, done = done)

    return (path_ins, outputs, options, spec)
    
### Here you should wait for the output. The output will comprise a file for each gene with the species sequence recovered.


##########################################################################################################################
###############################################---- MAFFT ----#############################################################
##########################################################################################################################

# Here go to folder 6.Retrieve and ls -1
# Get the gene names and write them in genes = []
# We found 3432 genes for Ceroxyloids using the Arecoideae target file.

def mafft(genes, path_in, path_out, done, gene):
    """Aligning all the sequences for each gene."""
    path_ins = [path_in+genes]
    outputs = [done, path_out+gene+"_aligned.fasta"] 
    options = {'cores': 4, 'memory': "4g", 'walltime': "4:00:00", 'account':"cryptocarya"}

    spec = """

    source /home/laurakf/miniconda3/etc/profile.d/conda.sh

    conda activate Mafft

    cd {path_in}

    mafft --thread 4 --globalpair --adjustdirectionaccurately --maxiterate 1000 {genes} > {path_out}{gene}_aligned.fasta

    touch {done}

    """.format(genes = genes, done = done, path_in = path_in, path_out = path_out, gene = gene)

    return (path_ins, outputs, options, spec)
    
#It is a good idea to rename the fasta files here.


########################################################################################################################
###############################################---- TRIMAL ----#########################################################
########################################################################################################################

#Cleaning according trimal
#Get raw alignments and trim them according to a gap threshold.

def gt_trimming(path_in, path_out, done, gene):
    """ Use trimal for trimming all alignments for each of the GT values specified"""
    inputs = [path_in+gene+"_aligned.fasta"]
    outputs = [done]
    options = {'cores': 1, 'memory': "5g", 'walltime': "1:00:00", 'account':"cryptocarya"}

    spec="""
    
    #Activating Trimal
    source /home/laurakf/miniconda3/etc/profile.d/conda.sh
    conda activate Trimal

    #Go to alignments folder
    cd {path_in}

    #Running gaptrimming.sh
    cd {path_out}
    bash /home/laurakf/cryptocarya/Scripts/gap_trimming.sh -g {gene}_aligned.fasta.old

    touch {done}

    """.format(path_in = path_in, done = done, gene = gene, path_out = path_out)

    return(inputs, outputs, options, spec)
    
    
########################################################################################################################
#################################################---- AMAS ----#########################################################
########################################################################################################################


######## Calculating amas summary (raw)

#For raw alignments
def amas_raw(path_in, done):
    """Creating summary files for all the trimmed alignments for each raw alignment"""
    inputs = [path_in]
    outputs = [path_in+"summary_0.txt", done]
    options = {'cores': 1, 'memory': "2g", 'walltime': "1:00:00", 'account':"cryptocarya"}

    spec="""


    #Activating AMAS
    source /home/laurakf/miniconda3/etc/profile.d/conda.sh
    conda activate Amas
    
    cd {path_in}

    #Calculating amas summary
    /home/laurakf/cryptocarya/Scripts/AMAS/amas/AMAS.py summary -f fasta -d dna -i *.fasta

    mv summary.txt summary_0.txt

    touch {done}
    
    """.format(path_in = path_in, done = done)

    return(inputs, outputs, options, spec)


# Hereafter you need to do some manual work and remove the headlines statistics.


######## Calculating amas summary (gt)

#For cutoff (gt) alignments
def amas_gt(path_in, cut_off, done):
    """Creating summary files for all the trimmed alignments for each raw alignment"""
    inputs = [path_in, path_in+cut_off]
    outputs = [path_in+"summary_0.1.txt", path_in+"summary_0.90.txt", done]
    options = {'cores': 1, 'memory': "2g", 'walltime': "1:00:00", 'account':"cryptocarya"}

    spec="""

    #Activating AMAS
    source /home/laurakf/miniconda3/etc/profile.d/conda.sh
    conda activate Amas
    
    cd {path_in}{cut_off} 

    #Calculating amas summary
    /home/laurakf/cryptocarya/Scripts/AMAS/amas/AMAS.py summary -f fasta -d dna -i *.fasta.old
   
    mv summary.txt summary_{cut_off}.txt 

    touch {done}
    
    """.format(path_in = path_in, cut_off = cut_off, done = done)

    return(inputs, outputs, options, spec)


# ########################################################################################################################
# #############################################---- Optrimal ----#########################################################
# ########################################################################################################################

# #Getting the best alignment for each gene 
# def optrim(path_in, path_out, done):
#     """Select the best alignments according to the gt value"""
#     inputs = [path_in+"0.1",path_in+"0.15",path_in+"0.20",path_in+"0.25",path_in+"0.30",
#     path_in+"0.35",path_in+"0.40",path_in+"0.45",path_in+"0.50",path_in+"0.55",path_in+"0.60",path_in+"0.65",
#     path_in+"0.70",path_in+"0.75",path_in+"0.80",path_in+"0.85",path_in+"0.90"]
#     outputs = [done, path_out+"optimal_final_results/"]
#     options = {'cores': 4, 'memory': "5g", 'walltime': "01:00:00", 'account':"cryptocarya"}

#     spec="""

#     #Going to folder with trimmed files
#     cd {path_in}

#     Rscript --vanilla /home/laurakf/cryptocarya/Scripts/optrimal.R

#     mv dldp_* {path_out}+"optrim_output/"
    
#     mv optimal_final_results {path_out}+"optimal_final_results/"

#     touch {done}

#     """.format(path_in = path_in, path_out = path_out, done = done)

#     return(inputs, outputs, options, spec)

# #Here we should include remove empty files?

# #HERE you can have an annoying error in R: "Error in real_loss[1:(length(real_loss) - 1)] :  only 0's may be mixed with negative subscripts"
# # Download all the summary and cutoff files.
# # Run the R code in your computer
# # Check the file lost, and see if there is a gene that only has 1.00000 for all gt values, and delete those
# # Copy the genes that you delete from raw_alignments to the trimal_output


# ##########################################################################################################################
# ##############################################---- CIALIGN ----###########################################################
# ##########################################################################################################################

# def cialign1(genes, path_in, path_out, done):
#     """Cleaning alignments using cialign default."""
#     inputs = [path_in + genes + "_aligned.fasta.old"]
#     outputs = [path_out+genes+"_cialign.fasta_cleaned.fasta", done]
#     options = {'cores': 8, 'memory': "100g", 'walltime': "12:00:00", 'account':"cryptocarya"}

#     spec = """
#     source activate CIAlign

#     cd {path_in}

#     CIAlign --infile {genes}_aligned.fasta.old --all --outfile_stem {path_out}{genes}_cialign.fasta

#     touch {done}

#     """.format(genes = genes, done = done, path_in = path_in, path_out = path_out)

#     return (inputs, outputs, options, spec)
    
# #Here we lost the gene 2683. This gene is too divergent

# ########################################################################################################################
# ###############################################---- TAPER ----##########################################################
# ########################################################################################################################

# def taper(path_in, genes, path_out, done):
#     """Using TAPER AFTER CIAlign to remove errors in small species-specific stretches of the multiple sequence alignments"""
#     inputs = [path_in+genes+"_cialign.fasta_cleaned.fasta", done]
#     outputs = ["/home/laurakf/cryptocarya/Workflow/Test/13_Taper/"+genes+"_output_tapper.fasta", done]
#     options = {'cores': 1, 'memory': "40g", 'walltime': "02:00:00", 'account':"cryptocarya"}

#     spec = """
     
#     cd {path_in}
        
#     #Activate the enviroment
#     source activate Taper
        
#     julia /home/laurakf/cryptocarya/Programs/TAPER-master/correction_multi.jl {genes}_cialign.fasta_cleaned.fasta > {genes}_output_tapper.fasta 
    
#     mv *_output_tapper.fasta {path_out}

#     touch {done}
        
#     """.format(path_in = path_in, genes = genes, path_out = path_out, done = done)

#     return (inputs, outputs, options, spec)
    

# ########################################################################################################################
# ##############################################---- IQTREE ----##########################################################
# ########################################################################################################################

# def iqtree(path_in, genes, done):
#     """Using IQTREE to construct a phylogenetic hypotheses for each gene"""
#     inputs = [path_in+genes+"_output_tapper.fasta", done]
#     outputs = [path_out+genes+"_output_tapper.fasta.treefile", done]
#     options = {'cores': 2, 'memory': "40g", 'walltime': "12:00:00", 'account':"cryptocarya"}

#     spec = """
     
#     cd {path_in}
        
#     #Activate the enviroment
#     source activate IQtree
        
#     iqtree2 -s {genes}_output_tapper.fasta -T AUTO -m MFP -B 1000 
    
#     mv *treefile {path_out}
#     mv *_output_tapper.fasta.model.gz {path_out}
#     mv *output_tapper.fasta.contree {path_out}
#     mv *output_tapper.fasta.bionj {path_out}
#     mv *output_tapper.fasta.ckp.gz {path_out}
#     mv *_output_tapper.fasta.iqtree {path_out}
#     mv *_output_tapper.fasta.log {path_out}
#     mv *tapper.fasta.mldist {path_out}
#     mv *tapper.fasta.splits.nex {path_out}
#     mv *output_tapper.fasta.uniqueseq.phy {path_out}

#     touch {done}
    
#     """.format(path_in = path_in, genes = genes, done = done)

#     return (inputs, outputs, options, spec) 
 
# # We also lost the gene: because it has less 4 species only, and it does not make sense to perform a bootstrap (Iqtree)
# #Iqtree_32
# #Iqtree_320
# #Iqtree_509
# #Iqtree_536
# #Iqtree_709
# #Iqtree_774
# #Iqtree_870
# #Iqtree_873
# #Iqtree_874
# #Iqtree_881
# #Iqtree_1111
# #Iqtree_1204
# #Iqtree_1421
# #Iqtree_1576
# #Iqtree_1654
# #Iqtree_1687
# #Iqtree_2480
# #Iqtree_2584
# #Iqtree_2661
# #Iqtree_2680
# #Iqtree_2973

# # Therefore we have removed all the genes cited + 2683 (too divergent and excluded by CIalign) 

# #Here you should stop and go to folder /home/paola/faststorage/17.Final_organization/5.Ceroxyloids/12.IQtree and do cat *treefile > gene_trees.nex

# ########################################################################################################################
# #####################################---- Astral Tree Search ----#####################################################
# ########################################################################################################################


# def astral_tapper(path_in, gene_tree_file, output, done):
#     """Using Astral to construct a species tree based on the genetrees"""
#     inputs = [path_in+"gene_trees.nex", done]
#     outputs = [path_in + output, done]
#     options = {'cores': 20, 'memory': "40g", 'walltime': "48:00:00", 'account':"cryptocarya"}

#     spec = """
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
                     
#     cd {path_in} 
  
#     java -jar /home/laurakf/cryptocarya/Programs/Astral/astral.5.7.8.jar -i {gene_tree_file} -o {output} 2> log_out.log
    
#     bash rename.sh

#     """.format(path_in = path_in, gene_tree_file = gene_tree_file, output = output, done = done)

#     return (inputs, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################

#Species removed from pipeline as they had no gene recovery.
sp = ["Ocotea-foetens-WE521","Ocotea-gabonensis-WE522","Ocotea-meziana-WE523","Pleurothyrium-cuneifolium-WE524","Mespilodaphne-cymbarum-WE525","Damburneya-gentlei-WE526","Ocotea-glaucosericea-WE527","Ocotea-complicata-WE528","Ocotea-javitensis-WE529","Ocotea-skutchii-WE530","Ocotea-sinuata-WE531","Ocotea-botrantha-WE532","Nectandra-lineatifolia-WE533"] 


for i in range(len(sp)):
    #### Running fastqc on raw data
    gwf.target_from_template('fastqc_raw_'+str(i), fastqc_raw(name = sp[i],
                                                        path_in= "/home/laurakf/cryptocarya/RawData/Test/",
                                                        path_out = "/home/laurakf/cryptocarya/Workflow/Test/01_FastQC/",
                                                        done = "/home/laurakf/cryptocarya/Workflow/Test/01_FastQC/done/"+sp[i]))

#### Running multiqc on raw data
gwf.target_from_template('multiqc_raw', multiqc_raw(path_in= "/home/laurakf/cryptocarya/Workflow/Test/01_FastQC/",
                                                    path_out = "/home/laurakf/cryptocarya/Workflow/Test/02_MultiQC/",
                                                    done = "/home/laurakf/cryptocarya/Workflow/Test/02_MultiQC/done/multiqc_raw"))

sp = ["Ocotea-foetens-WE521","Ocotea-gabonensis-WE522","Ocotea-meziana-WE523","Pleurothyrium-cuneifolium-WE524","Mespilodaphne-cymbarum-WE525","Damburneya-gentlei-WE526","Ocotea-glaucosericea-WE527","Ocotea-complicata-WE528","Ocotea-javitensis-WE529","Ocotea-skutchii-WE530","Ocotea-sinuata-WE531","Ocotea-botrantha-WE532","Nectandra-lineatifolia-WE533"] 

for i in range(len(sp)):
    #### Running Trimmomatic
    gwf.target_from_template('trimmomatic_'+str(i), trimmomatic(name = sp[i],
                                                        path_in= "/home/laurakf/cryptocarya/Workflow/Test/01_FastQC/data/", 
                                                        path_out = "/home/laurakf/cryptocarya/Workflow/Test/03_Trimmomatic/slidingwindow/",
                                                        done = "/home/laurakf/cryptocarya/Workflow/Test/03_Trimmomatic/slidingwindow/done/"+sp[i]))

    # #### Running fastqc on the trimmed data
    gwf.target_from_template('fastqc_trimmed_'+str(i), fastqc_trimmed(name = sp[i],
                                                        path_in= "/home/laurakf/cryptocarya/Workflow/Test/03_Trimmomatic/slidingwindow/secapr_postrim/",
                                                        path_out = "/home/laurakf/cryptocarya/Workflow/Test/04_FastQC/slidingwindow/",
                                                        done = "/home/laurakf/cryptocarya/Workflow/Test/04_FastQC/slidingwindow/done/"+sp[i]))                                                   


# #### Running multiqc on trimmed data
# gwf.target_from_template('multiqc_trimmed', multiqc_trimmed(path_in= "/home/laurakf/cryptocarya/Workflow/Test/04_FastQC/slidingwindow/",
#                                                     path_out = "/home/laurakf/cryptocarya/Workflow/Test/05_MultiQC/slidingwindow/",
#                                                     done = "/home/laurakf/cryptocarya/Workflow/Test/05_MultiQC/slidingwindow/done/multiqc_trimmed"))


sp = ["Ocotea-foetens-WE521","Ocotea-gabonensis-WE522","Ocotea-meziana-WE523","Pleurothyrium-cuneifolium-WE524","Mespilodaphne-cymbarum-WE525","Damburneya-gentlei-WE526","Ocotea-glaucosericea-WE527","Ocotea-complicata-WE528","Ocotea-javitensis-WE529","Ocotea-skutchii-WE530","Ocotea-sinuata-WE531"] 
# Taken "Ocotea-botrantha-WE532" and "Nectandra-lineatifolia-WE533" out. They do not seem to work. 

for i in range(len(sp)):
    
    #### Running Hybpiper
    gwf.target_from_template('Hybpiper_'+str(i), hybpiper(name = sp[i],
                                                        p1 = "_1P.fastq",
                                                        p2 = "_2P.fastq",
                                                        un = "_UN.fastq",
                                                        path_out= "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/",
                                                        path_in = "/home/laurakf/cryptocarya/Workflow/Test/03_Trimmomatic/slidingwindow/",
                                                        done = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/HybPiper/"+sp[i]))

#### Getting stats and heatmap
gwf.target_from_template('stats', stats(path_out= "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/Stats_Heatmap/",
                                                path_in = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/",
                                                in_done = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/HybPiper/"+sp[i],
                                                name = sp[i],
                                                done = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/Stats/"+sp[i]))
                                                  
    # #### Paralogs
    # gwf.target_from_template('Paralogs_'+str(i), paralogs(name = sp[i],
    #                                                     path_in = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/",
    #                                                     done = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/Paralogs/"+sp[i],
    #                                                     in_done="/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/HybPiper/"+sp[i]))


sp = ["Ocotea-foetens-WE521","Ocotea-gabonensis-WE522","Ocotea-meziana-WE523","Pleurothyrium-cuneifolium-WE524","Mespilodaphne-cymbarum-WE525","Damburneya-gentlei-WE526","Ocotea-glaucosericea-WE527","Ocotea-complicata-WE528","Ocotea-javitensis-WE529","Ocotea-skutchii-WE530","Ocotea-sinuata-WE531"] 
# Taken "Ocotea-botrantha-WE532" and "Nectandra-lineatifolia-WE533" out. They do not seem to work. 

## paralogs
for i in range(len(sp)):
    #### Paralogs
    if os.path.isfile("/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/"+sp[i]+"/"+sp[i]+"_genes_with_long_paralog_warnings.txt"):
        gwf.target_from_template('Paralogs_'+str(i), paralogs(name = sp[i],
                                                            path_out = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/Paralogs/",
                                                            path_in = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/",
                                                            in_done = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/HybPiper/"+sp[i],
                                                            done = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/Paralogs/"+sp[i]))

    ## No paralogs
    else:
        gwf.target_from_template('No_Paralogs_'+str(i), no_paralogs(name = sp[i],
                                                                path_in = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/",
                                                                done = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/Paralogs/"+sp[i],
                                                                in_done = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/HybPiper/"+sp[i],
                                                                no_paralogs="/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/No_Paralogs/"+sp[i]))                   

    #### Coverage
    gwf.target_from_template('Coverage_'+str(i), coverage(name = sp[i],
                                                        path_in = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/",
                                                        all_bam = "_all.bam",
                                                        all_sorted_bam ="_all_sorted.bam",
                                                        all_sorted_bam_bai="_all_sorted.bam.bai",
                                                        bam =".bam",
                                                        cov=".cov",
                                                        fasta = ".fasta",
                                                        fasta_amb = ".fasta.amb",
                                                        fasta_ann = ".fasta.ann",
                                                        fasta_bwt = ".fasta.bwt",
                                                        fasta_pac = ".fasta.pac",
                                                        fasta_sa = ".fasta.sa",
                                                        trimmed_fasta = "_trimmed.fasta",
                                                        up_bam = "_up.bam",
                                                        path_out = "/home/laurakf/cryptocarya/Workflow/Test/07_Coverage/",
                                                        done = "/home/laurakf/cryptocarya/Workflow/Test/07_Coverage/done/Coverage/"+sp[i],
                                                        dir_wrk = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/",
                                                        dir_in ="/home/laurakf/cryptocarya/Workflow/Test/03_Trimmomatic/slidingwindow/", #Folder with clean reads + unpaired
                                                        dir_out ="/home/laurakf/cryptocarya/Workflow/Test/07_Coverage/")) # folder with coverage

#### Retrieve sequences and sort into files with gene names
gwf.target_from_template('retrieve', retrieve(path_in ="/home/laurakf/cryptocarya/Workflow/Test/07_Coverage/", 
                                              done = "/home/laurakf/cryptocarya/Workflow/Test/08_Retrieve/done/retrieve"))


gene = ["4471", "4527", "4691", "4724", "4744", "4757", "4793", "4796", "4802", "4806", "4848", "4889", "4890", "4893", "4932", "4942", "4951", "4954", "4989", "4992", "5018", "5032", "5034", "5038", "5090", "5104", "5116", "5123", "5131", "5138", "5162", "5163", "5168", "5177", "5188", "5200", "5206", "5220", "5257", "5264", "5271", "5273", "5280", "5296", "5299", "5304", "5318", "5326", "5328", "5333", "5335", "5339", "5343", "5348", "5354", "5355", "5366", "5398", "5404", "5406", "5421", "5426", "5427", "5428", "5430", "5434", "5449", "5454", "5460", "5463", "5464", "5469", "5477", "5489", "5502", "5513", "5528", "5531", "5536", "5551", "5554", "5562", "5578", "5594", "5596", "5599", "5614", "5620", "5634", "5639", "5644", "5656", "5660", "5664", "5670", "5699", "5702", "5703", "5716", "5721", "5744", "5770", "5772", "5791", "5802", "5815", "5816", "5821", "5822", "5840", "5841", "5842", "5843", "5849", "5853", "5857", "5858", "5859", "5865", "5866", "5870", "5893", "5894", "5899", "5910", "5913", "5918", "5919", "5921", "5926", "5933", "5940", "5942", "5943", "5944", "5945", "5949", "5950", "5958", "5960", "5968", "5974", "5977", "5980", "5981", "5990", "6000", "6003", "6004", "6016", "6026", "6029", "6034", "6036", "6038", "6041", "6048", "6050", "6051", "6056", "6064", "6068", "6072", "6098", "6110", "6114", "6119", "6128", "6130", "6139", "6148", "6164", "6175", "6176", "6198", "6216", "6221", "6226", "6227", "6238", "6258", "6265", "6274", "6282", "6284", "6295", "6298", "6299", "6303", "6318", "6320", "6363", "6366", "6373", "6376", "6378", "6383", "6384", "6387", "6389", "6393", "6398", "6404", "6405", "6406", "6407", "6412", "6420", "6432", "6439", "6447", "6450", "6454", "6457", "6458", "6459", "6460", "6462", "6483", "6487", "6488", "6492", "6494", "6496", "6498", "6500", "6506", "6507", "6526", "6527", "6528", "6532", "6533", "6538", "6540", "6544", "6550", "6552", "6559", "6563", "6570", "6572", "6601", "6620", "6631", "6636", "6639", "6641", "6649", "6652", "6660", "6667", "6685", "6689", "6713", "6717", "6732", "6733", "6738", "6746", "6779", "6782", "6785", "6792", "6797", "6825", "6848", "6854", "6859", "6860", "6865", "6875", "6882", "6883", "6909", "6913", "6914", "6924", "6933", "6946", "6947", "6954", "6955", "6958", "6961", "6962", "6968", "6978", "6979", "6992", "6995", "7021", "7024", "7029", "7067", "7111", "7128", "7135", "7136", "7141", "7174", "7194", "7241", "7273", "7279", "7313", "7324", "7325", "7331", "7333", "7336", "7363", "7367", "7371", "7572", "7577", "7583", "7602", "7628"]
genes = ["4471.FNA", "4527.FNA", "4691.FNA", "4724.FNA", "4744.FNA", "4757.FNA", "4793.FNA", "4796.FNA", "4802.FNA", "4806.FNA", "4848.FNA", "4889.FNA", "4890.FNA", "4893.FNA", "4932.FNA", "4942.FNA", "4951.FNA", "4954.FNA", "4989.FNA", "4992.FNA", "5018.FNA", "5032.FNA", "5034.FNA", "5038.FNA", "5090.FNA", "5104.FNA", "5116.FNA", "5123.FNA", "5131.FNA", "5138.FNA", "5162.FNA", "5163.FNA", "5168.FNA", "5177.FNA", "5188.FNA", "5200.FNA", "5206.FNA", "5220.FNA", "5257.FNA", "5264.FNA", "5271.FNA", "5273.FNA", "5280.FNA", "5296.FNA", "5299.FNA", "5304.FNA", "5318.FNA", "5326.FNA", "5328.FNA", "5333.FNA", "5335.FNA", "5339.FNA", "5343.FNA", "5348.FNA", "5354.FNA", "5355.FNA", "5366.FNA", "5398.FNA", "5404.FNA", "5406.FNA", "5421.FNA", "5426.FNA", "5427.FNA", "5428.FNA", "5430.FNA", "5434.FNA", "5449.FNA", "5454.FNA", "5460.FNA", "5463.FNA", "5464.FNA", "5469.FNA", "5477.FNA", "5489.FNA", "5502.FNA", "5513.FNA", "5528.FNA", "5531.FNA", "5536.FNA", "5551.FNA", "5554.FNA", "5562.FNA", "5578.FNA", "5594.FNA", "5596.FNA", "5599.FNA", "5614.FNA", "5620.FNA", "5634.FNA", "5639.FNA", "5644.FNA", "5656.FNA", "5660.FNA", "5664.FNA", "5670.FNA", "5699.FNA", "5702.FNA", "5703.FNA", "5716.FNA", "5721.FNA", "5744.FNA", "5770.FNA", "5772.FNA", "5791.FNA", "5802.FNA", "5815.FNA", "5816.FNA", "5821.FNA", "5822.FNA", "5840.FNA", "5841.FNA", "5842.FNA", "5843.FNA", "5849.FNA", "5853.FNA", "5857.FNA", "5858.FNA", "5859.FNA", "5865.FNA", "5866.FNA", "5870.FNA", "5893.FNA", "5894.FNA", "5899.FNA", "5910.FNA", "5913.FNA", "5918.FNA", "5919.FNA", "5921.FNA", "5926.FNA", "5933.FNA", "5940.FNA", "5942.FNA", "5943.FNA", "5944.FNA", "5945.FNA", "5949.FNA", "5950.FNA", "5958.FNA", "5960.FNA", "5968.FNA", "5974.FNA", "5977.FNA", "5980.FNA", "5981.FNA", "5990.FNA", "6000.FNA", "6003.FNA", "6004.FNA", "6016.FNA", "6026.FNA", "6029.FNA", "6034.FNA", "6036.FNA", "6038.FNA", "6041.FNA", "6048.FNA", "6050.FNA", "6051.FNA", "6056.FNA", "6064.FNA", "6068.FNA", "6072.FNA", "6098.FNA", "6110.FNA", "6114.FNA", "6119.FNA", "6128.FNA", "6130.FNA", "6139.FNA", "6148.FNA", "6164.FNA", "6175.FNA", "6176.FNA", "6198.FNA", "6216.FNA", "6221.FNA", "6226.FNA", "6227.FNA", "6238.FNA", "6258.FNA", "6265.FNA", "6274.FNA", "6282.FNA", "6284.FNA", "6295.FNA", "6298.FNA", "6299.FNA", "6303.FNA", "6318.FNA", "6320.FNA", "6363.FNA", "6366.FNA", "6373.FNA", "6376.FNA", "6378.FNA", "6383.FNA", "6384.FNA", "6387.FNA", "6389.FNA", "6393.FNA", "6398.FNA", "6404.FNA", "6405.FNA", "6406.FNA", "6407.FNA", "6412.FNA", "6420.FNA", "6432.FNA", "6439.FNA", "6447.FNA", "6450.FNA", "6454.FNA", "6457.FNA", "6458.FNA", "6459.FNA", "6460.FNA", "6462.FNA", "6483.FNA", "6487.FNA", "6488.FNA", "6492.FNA", "6494.FNA", "6496.FNA", "6498.FNA", "6500.FNA", "6506.FNA", "6507.FNA", "6526.FNA", "6527.FNA", "6528.FNA", "6532.FNA", "6533.FNA", "6538.FNA", "6540.FNA", "6544.FNA", "6550.FNA", "6552.FNA", "6559.FNA", "6563.FNA", "6570.FNA", "6572.FNA", "6601.FNA", "6620.FNA", "6631.FNA", "6636.FNA", "6639.FNA", "6641.FNA", "6649.FNA", "6652.FNA", "6660.FNA", "6667.FNA", "6685.FNA", "6689.FNA", "6713.FNA", "6717.FNA", "6732.FNA", "6733.FNA", "6738.FNA", "6746.FNA", "6779.FNA", "6782.FNA", "6785.FNA", "6792.FNA", "6797.FNA", "6825.FNA", "6848.FNA", "6854.FNA", "6859.FNA", "6860.FNA", "6865.FNA", "6875.FNA", "6882.FNA", "6883.FNA", "6909.FNA", "6913.FNA", "6914.FNA", "6924.FNA", "6933.FNA", "6946.FNA", "6947.FNA", "6954.FNA", "6955.FNA", "6958.FNA", "6961.FNA", "6962.FNA", "6968.FNA", "6978.FNA", "6979.FNA", "6992.FNA", "6995.FNA", "7021.FNA", "7024.FNA", "7029.FNA", "7067.FNA", "7111.FNA", "7128.FNA", "7135.FNA", "7136.FNA", "7141.FNA", "7174.FNA", "7194.FNA", "7241.FNA", "7273.FNA", "7279.FNA", "7313.FNA", "7324.FNA", "7325.FNA", "7331.FNA", "7333.FNA", "7336.FNA", "7363.FNA", "7367.FNA", "7371.FNA", "7572.FNA", "7577.FNA", "7583.FNA", "7602.FNA", "7628.FNA"]


#### MAFFT
for i in range(len(genes)):
    gwf.target_from_template('Mafft_'+str(i), mafft(genes = genes[i],
                                                        path_out= "/home/laurakf/cryptocarya/Workflow/Test/09_Mafft/",
                                                        path_in = "/home/laurakf/cryptocarya/Workflow/Test/08_Retrieve/",
                                                        gene = gene[i],
                                                        done = "/home/laurakf/cryptocarya/Workflow/Test/09_Mafft/done/"+genes[i]))


#### Trimal according to a pre-defined gt values
for i in range(len(gene)):
   gwf.target_from_template('gt_trimming_'+gene[i], gt_trimming(gene = gene[i],
                                                        path_in = "/home/laurakf/cryptocarya/Workflow/Test/09_Mafft/",
                                                        path_out = "/home/laurakf/cryptocarya/Workflow/Test/10_Trimal/",
                                                        done = "/home/laurakf/cryptocarya/Workflow/Test/10_Trimal/done/"+gene[i]))

#### Generating AMAS statistics for raw_alignments
gwf.target_from_template('amas_raw', amas_raw(path_in = "/home/laurakf/cryptocarya/Workflow/Test/10_Trimal/",
                                        done = "/home/laurakf/cryptocarya/Workflow/Test/10_Trimal/done/AMAS_raw/raw"))

cut_off = ["0.1", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45", "0.50", "0.55", "0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90"]

#### Generating AMAS statistics for gt_alignments
for i in range(len(cut_off)):
    gwf.target_from_template('amas_gt_'+cut_off[i], amas_gt(path_in = "/home/laurakf/cryptocarya/Workflow/Test/10_Trimal/",
                                                cut_off = cut_off[i],
                                                done = "/home/laurakf/cryptocarya/Workflow/Test/10_Trimal/done/AMAS_gt"+cut_off[i]))


# #### Optrimal
# gwf.target_from_template('optrim', optrim(path_in = "/home/laurakf/cryptocarya/Workflow/Test/10_Trimal/",
#                                                          done = "optimal_final_results /home/laurakf/cryptocarya/Workflow/Test/11_Optrimal/done/"+genes[i], 
#                                                          path_out = "/home/laurakf/cryptocarya/Workflow/Test/11_Optrimal/"))

                                               
# # Running CIAlign on the trimmed_fasta - Including Paralogs
# for i in range(0, len(genes)):
#     gwf.target_from_template('Cialign1'+str(i), cialign1(genes = genes[i],
#                                               path_in = ""/home/laurakf/cryptocarya/Workflow/Test/11_Optrimal/optimal_final_results/",
#                                               path_out = "/home/laurakf/cryptocarya/Workflow/Test/12_CIAlign/",
#                                               done = "/home/laurakf/cryptocarya/Workflow/Test/12_CIAlign/done/"+genes[i]))

# ## Running TAPER after CIALIGN
# for i in range(0, len(genes)):
#     gwf.target_from_template('Taper_'+str(i), taper(genes = genes[i],
#                                                     path_in = "/home/laurakf/cryptocarya/Workflow/Test/12_CIAlign/",
#                                                     path_out = "/home/laurakf/cryptocarya/Workflow/Test/13_Taper/",
#                                                     done = "/home/laurakf/cryptocarya/Workflow/Test/13_Taper/done"+genes[i]))
                                                    
# #Running IQTREE for files trimmed with trimal and CIAlign                                             
# for i in range(0, len(genes)):
#    gwf.target_from_template('Iqtree_'+str(i), iqtree(genes = genes[i],
#                                                     path_in = "/home/laurakf/cryptocarya/Workflow/Test/13_Taper/"))  
                                                    
# # Running ASTRAL 
# gwf.target_from_template('astral_tapper', astral_tapper(path_in = "/home/paola/faststorage/17.Final_organization/5.Ceroxyloids/13.Astral/",
#                                                     gene_tree_file="gene_trees.nex",
#                                                     output="astral_tree_only_posterior_probability.tre"))



