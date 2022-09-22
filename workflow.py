
# ------------------------------------------------------------------------------------------------------------------------
# This workflow is used to transform the raw sequence data into sequences ready for alignment.
# Workflow is the following

# 1: Secapr quality check of the sequences
# 2: Trimming of the sequences using trimmomatic
# 3: Hybpipering the species in order to create the exons of the baits for each species
# 4: Checking for Paralogs
# 5: Running Intronerate again to get the Introns for each species
# 6: Trimming for Coverage of sequencing and joining the exon and intron of each species to create supercontigs

# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# Author: Oscar Wrisberg
# Date: 10/11/2021
# ------------------------------------------------------------------------------------------------------------------------

from os import O_SYNC, name
from gwf import Workflow
import os.path
# import math
# import glob

gwf = Workflow()

########################################################################################################################
################################################---- Fastqc quality check raw ----#######################################
########################################################################################################################
def fastqc_raw(name,path_in ,path_out, done,):
    """Quality checking using fastqc as this should work on individual species"""
    inputs = []
    outputs = [path_out+name+"_R1_fastqc.html",path_out+name+"_R2_fastqc.html", done]
    options = {'cores': 1, 'memory': "8g", 'walltime': "00:30:00", 'account':"cryptocarya"}


    spec = """

    echo {name}

    /faststorage/project/cryptocarya/Programs/FastQC/fastqc -o {path_out} {path_in}{name}:_R1.fastq {path_in}{name}:_R2.fastq
    
    echo touching {done}
    touch {done}

    """.format(path_in = path_in,name = name, path_out = path_out, done = done)

    return (inputs, outputs, options, spec)

########################################################################################################################
################################################---- Fastqc quality check trimmed ----#######################################
########################################################################################################################
# def fastqc_trimmed(species,path_in ,path_out, done,):
#     """Quality checking using fastqc as this should work on individual species"""
#     inputs = [path_in+species+"_UN.fastq", path_in+species+"_1PU.fastq", path_in+species+"_2PU.fastq","/home/owrisberg/Coryphoideae/work_flow/02_trimmed/done/"+species]
#     outputs = [path_out+species+"_1PU_fastqc.html", path_out+species+"_2PU_fastqc.html",path_out+species+"_UN_fastqc.html" ,done]
#     options = {'cores': 1, 'memory': "10g", 'walltime': "00:30:00", 'account':"Coryphoideae"}


#     spec = """
#     source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
#     conda activate secapr_env

#     fastqc -o {output} {path_in}{species}_1PU.fastq {path_in}{species}_2PU.fastq {path_in}{species}_UN.fastq
    
#     touch {done}

#     """.format(path_in = path_in,species = species, output = path_out, done = done)

#     return (inputs, outputs, options, spec)


# # ########################################################################################################################
# # ################################################---- Trimmomatic ----###################################################
# # ########################################################################################################################
# def trimmomatic(species, path_in, path_out, done):
#     """Trimming raw data using trimmomatic with custom adapter.
#     Afterwards combines paired and unpaired reads for forward and reverse reads respectively for each species 
#     to enable post-trimming secapr quality_check for comparability before and after trimming """
#     inputs = []
#     outputs = [path_out+species+"_UN.fastq",path_out+"secapr_postrim/"+species+"_UN.fastq", done, path_out+species+"_1P.fastq", path_out+species+"_2P.fastq"]
#     options = {'cores': 16, 'memory': "10g", 'walltime': "01:00:00", 'account':"Coryphoideae"}

#     spec = """
#     source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
#     conda activate trimmomatic_env

#     trimmomatic PE -threads 16 -phred33 {input}_R1.fastq {input}_R2.fastq -baseout {output}.fastq\
#     ILLUMINACLIP:/home/owrisberg/miniconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10:1:true\
#     LEADING:3\
#     TRAILING:3\
#     MAXINFO:40:0.8\
#     MINLEN:36\
#     2>> stderr_trim_loop_output.txt

#     echo combining {path_out}{species}_1P.fastq and {path_out}{species}_1U.fastq into {path_out}secapr_postrim/{species}_1PU.fastq 
#     cat {path_out}{species}_1P.fastq {path_out}{species}_1U.fastq > {path_out}secapr_postrim/{species}_1PU.fastq 

#     echo combining {path_out}{species}_2P.fastq and {path_out}{species}_2U.fastq into {path_out}secapr_postrim/{species}_2PU.fastq 
#     cat {path_out}{species}_2P.fastq {path_out}{species}_2U.fastq > {path_out}secapr_postrim/{species}_2PU.fastq

#     echo combining {path_out}{species}_1U.fastq {path_out}{species}_2U.fastq > {path_out}{species}_UN.fastq
#     cat {path_out}{species}_1U.fastq {path_out}{species}_2U.fastq > {path_out}{species}_UN.fastq
#     cp {path_out}{species}_UN.fastq {path_out}secapr_postrim/


#     echo Removing {path_out}{species}_1U.fastq
#     rm {path_out}{species}_1U.fastq

#     echo Removing {path_out}{species}_2U.fastq
#     rm {path_out}{species}_2U.fastq

#     touch {done}

#     """.format(input = path_in + species, output = path_out+species, done = done, species = species, path_out = path_out)

#     return (inputs, outputs, options, spec)


# ########################################################################################################################
# ################################################---- Hybpiper ----######################################################
# ########################################################################################################################
# def hybpiper(species, p1, p2, un, path_out, path_in, done):
#     """Hybpiper."""
#     inputs = [path_in + species +p1, path_in + species + p2, path_in + species + un] # The files which the job will look for before it runs
#     outputs = [path_out + species, done] # The files which will have to be created in order for the job to be "completed"
#     options = {'cores': 1, 'memory': "20g", 'walltime': "100:00:00", 'account':"Coryphoideae"} #Slurm commands

#     spec = """

#     source /home/owrisberg/miniconda3/etc/profile.d/conda.sh

#     conda activate base

#     cd /scratch/$SLURM_JOBID
        
#     /home/owrisberg/Coryphoideae/github_code/HybPiper/reads_first.py --cpu 1 --readfiles {p1} {p2} --unpaired {un} -b /home/owrisberg/Coryphoideae/target_sequence/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta --prefix {species} --bwa

#     cp --recursive --update {species} /home/owrisberg/Coryphoideae/work_flow/03_hybpiper/

#     touch {done}
#     touch {out}{species}

#     """.format(species=species, p1 = path_in + species + p1,p2 = path_in + species + p2, un = path_in + species + un , out = path_out, done = done)


#     return (inputs, outputs, options, spec)

# ########################################################################################################################
# #############################################---- Paralogs ----#########################################################
# ########################################################################################################################

# def paralogs(species,path_in, done, no_paralogs, in_done):
#     """Find Paralog genes and write them in the file called paralog.txt"""
#     inputs = [path_in + species, in_done]
#     outputs = [done]
#     options = {'cores': 2, 'memory': "10g", 'walltime': "0:30:00", 'account':"Coryphoideae"}

#     spec = """
#     source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
#     conda activate base
    
#     if test -f /home/owrisberg/Coryphoideae/work_flow/03_hybpiper/{sp}/genes_with_paralog_warnings.txt; then
#         echo "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/{sp}/genes_with_paralog_warnings.txt exists" 
#         cd {path_in}
#         python /home/owrisberg/Coryphoideae/github_code/HybPiper/paralog_investigator.py {sp} 2>> paralog.txt
#     else
#         echo "the genes_with_paralog_warnings.txt does not exist and we run the no parallels part"
#         touch {np}
#     fi
    
#     touch {done}

#     """.format(sp = species, done = done, path_in = path_in, np = no_paralogs)
#     return (inputs, outputs, options, spec)

# def no_paralogs(species, path_in, done, no_paralogs):
#     """Wrapper script to continue pipeline when Hybpiper finds no paralogs"""
#     inputs = [path_in + species]
#     outputs = [done]
#     options = {'cores': 2, 'memory': "10g", 'walltime': "0:05:00", 'account':"Coryphoideae"}

#     spec = """

#     touch {done}
#     touch {np}

#     """.format(done=done, np=no_paralogs)
#     return(inputs, outputs, options, spec)

# # ########################################################################################################################
# # #############################################---- Intronerate ----######################################################
# # ########################################################################################################################

# def intronerate(species, path_in, done):
#     """Intronerate the sequencec from hybpiper."""
#     inputs = [path_in + species, path_in+"done/Hybpiper/"+species]
#     outputs = [done]
#     options = {'cores': 4, 'memory': "20g", 'walltime': "16:00:00", 'account':"Coryphoideae"}

#     spec = """
#     source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
#     conda activate base

#     cd {path_in}

#     python3 /home/owrisberg/Coryphoideae/github_code/HybPiper/intronerate.py --prefix {sp} &>> intronerate_out.txt
        
    
#     touch {done}

#     """.format(sp = species, done = done, path_in = path_in)

#     return (inputs, outputs, options, spec)

# # ########################################################################################################################
# # #############################################---- Coverage ----#########################################################
# # ########################################################################################################################
# def coverage(species, path_in, path_out, done,all_bam,all_sorted_bam, all_sorted_bam_bai, bam, cov,fasta,fasta_amb,fasta_ann,fasta_bwt,fasta_pac,fasta_sa,trimmed_fasta,up_bam,dir_in,dir_out):
#     """Calculating coverage of sequences."""
#     inputs = [path_in+species, path_in+"done/Intronerate/"+species]
#     outputs = [path_out+species+all_bam,
#      path_out+species+all_sorted_bam,
#       path_out+species+all_sorted_bam_bai,
#        path_out+species+bam,
#     path_out+species+cov,
#      path_out+species+fasta,
#       path_out+species+fasta_amb,
#        path_out+species+fasta_ann,
#         path_out+species+fasta_bwt,

#     path_out+species+fasta_pac,
#      path_out+species+fasta_sa,
#       path_out+species+trimmed_fasta,
#        path_out+species+up_bam,done] #ALL the output files
#     options = {'cores': 4, 'memory': "20g", 'walltime': "08:00:00", 'account':"Coryphoideae"}

#     spec = """
#     source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
#     conda activate base
    
#     cd {path_in}

#     python3 /home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/coverage.py {sp} {dir_in} {dir_out}
    
#     touch {done}

#     """.format(sp = species, done = done, path_in = path_in, dir_in = dir_in, dir_out = dir_out)

#     return (inputs, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################

sp = ["Ocotea-foetens-WE521","Ocotea-gabonensis-WE522","Ocotea-meziana-WE523","Pleurothyrium-cuneifolium-WE524","Mespilodaphne-cymbarum-WE525","Damburneya-gentlei-WE526","Ocotea-glaucosericea-WE527","Ocotea-complicata-WE528","Ocotea-javitensis-WE529","Ocotea-skutchii-WE530","Ocotea-sinuata-WE531","Ocotea-botrantha-WE532","Nectandra-lineatifolia-WE533"] 
#Species removed from pipeline as they had no gene recovery [3050,3188,3272,3300,3316, 3364, 3392, 3394]

for i in range(len(sp)):
    #### Running fastqc on raw data
    gwf.target_from_template('fastqc_raw_'+sp[i], fastqc_raw(name = sp[i],
                                                        path_in= "/home/laurakf/cryptocarya/RawData/Test/",
                                                        path_out = "/home/laurakf/cryptocarya/Workflow/Test/01_FastQC/",
                                                        done = "/home/laurakf/cryptocarya/Workflow/Test/01_FastQC/done"+sp[i]))


    # #### Running Trimmomatic
    # gwf.target_from_template('trimmomatic_'+sp[i], trimmomatic(species = sp[i],
    #                                                     path_in= "/home/owrisberg/Coryphoideae/work_flow/01_data/",
    #                                                     path_out = "/home/owrisberg/Coryphoideae/work_flow/02_trimmed/",
    #                                                     done = "/home/owrisberg/Coryphoideae/work_flow/02_trimmed/done/"+sp[i]))

    # #### Running fastqc on the trimmed data
    # gwf.target_from_template('fastqc_trimmed_'+sp[i], fastqc_trimmed(species = sp[i],
    #                                                     path_in= "/home/owrisberg/Coryphoideae/work_flow/02_trimmed/secapr_postrim/",
    #                                                     path_out = "/home/owrisberg/Coryphoideae/work_flow/00_secapr/1_trimmed/",
    #                                                     done = "/home/owrisberg/Coryphoideae/work_flow/00_secapr/done/trimmed_data/"+sp[i]))                                                   

    # #### Running Hybpiper
    # gwf.target_from_template('Hybpiper_'+sp[i], hybpiper(species = sp[i],
    #                                                     p1 = "_1P.fastq",
    #                                                     p2 = "_2P.fastq",
    #                                                     un = "_UN.fastq",
    #                                                     path_out= "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
    #                                                     path_in = "/home/owrisberg/Coryphoideae/work_flow/02_trimmed/",
    #                                                     done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Hybpiper/"+sp[i],))
                                                                      

    # #### Paralogs
    
    # gwf.target_from_template('Paralogs_'+sp[i], paralogs(species = sp[i],
    #                                                     path_in = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
    #                                                     done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Paralogs/"+sp[i],
    #                                                     no_paralogs="/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/No_paralogs/"+sp[i],
    #                                                     in_done="/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Hybpiper/"+sp[i]))
    # # else:
    # #     gwf.target_from_template('No_Paralogs_'+sp[i], no_paralogs(species = sp[i],
    # #                                                             path_in = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
    # #                                                             done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Paralogs/"+sp[i],
    # #                                                             no_paralogs="/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/No_paralogs/"+sp[i]))
     
    
    # #### Getting introns
    # gwf.target_from_template('Intronerate_'+sp[i], intronerate(species= sp[i],
    #                                                     path_in = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
    #                                                     done = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Intronerate/"+sp[i]))


    # #### Coverage
    # gwf.target_from_template('Coverage_'+sp[i], coverage(species = sp[i],
    #                                                     path_in = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
    #                                                     all_bam = "_all.bam",
    #                                                     all_sorted_bam ="_all_sorted.bam",
    #                                                     all_sorted_bam_bai="_all_sorted.bam.bai",
    #                                                     bam =".bam",
    #                                                     cov=".cov",
    #                                                     fasta = ".fasta",
    #                                                     fasta_amb = ".fasta.amb",
    #                                                     fasta_ann = ".fasta.ann",
    #                                                     fasta_bwt = ".fasta.bwt",
    #                                                     fasta_pac = ".fasta.pac",
    #                                                     fasta_sa = ".fasta.sa",
    #                                                     trimmed_fasta = "_trimmed.fasta",
    #                                                     up_bam = "_up.bam",
    #                                                     path_out = "/home/owrisberg/Coryphoideae/work_flow/04_coverage/",
    #                                                     done = "/home/owrisberg/Coryphoideae/work_flow/04_coverage/done/Coverage/"+sp[i],
    #                                                     dir_in ="/home/owrisberg/Coryphoideae/work_flow/02_trimmed/", #Folder with raw reads
    #                                                     dir_out ="/home/owrisberg/Coryphoideae/work_flow/04_coverage/")) # folder with coverage

