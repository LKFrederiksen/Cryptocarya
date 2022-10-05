
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
# Author: Laura Kragh Frederiksen (adapted from Oscar Wrisberg)
# Date: 
# ------------------------------------------------------------------------------------------------------------------------

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
def multiqc_trimmed(path_in ,path_out, done,):
    """Quality checking using multiqc"""
    inputs = [path_in+"Ocotea-foetens-WE521_1PU_fastqc.html",path_in+"Ocotea-gabonensis-WE522_1PU_fastqc.html",path_in+"Ocotea-meziana-WE523_1PU_fastqc.html",path_in+"Pleurothyrium-cuneifolium-WE524_1PU_fastqc.html",path_in+"Mespilodaphne-cymbarum-WE525_1PU_fastqc.html",path_in+"Damburneya-gentlei-WE526_1PU_fastqc.html",path_in+"Ocotea-glaucosericea-WE527_1PU_fastqc.html",path_in+"Ocotea-complicata-WE528_1PU_fastqc.html",path_in+"Ocotea-javitensis-WE529_1PU_fastqc.html",path_in+"Ocotea-skutchii-WE530_1PU_fastqc.html",path_in+"Ocotea-sinuata-WE531_1PU_fastqc.html",path_in+"Ocotea-botrantha-WE532_1PU_fastqc.html",path_in+"Nectandra-lineatifolia-WE533_1PU_fastqc.html",path_in+"Ocotea-foetens-WE521_2PU_fastqc.html",path_in+"Ocotea-gabonensis-WE522_2PU_fastqc.html",path_in+"Ocotea-meziana-WE523_2PU_fastqc.html",path_in+"Pleurothyrium-cuneifolium-WE524_2PU_fastqc.html",path_in+"Mespilodaphne-cymbarum-WE525_2PU_fastqc.html",path_in+"Damburneya-gentlei-WE526_2PU_fastqc.html",path_in+"Ocotea-glaucosericea-WE527_2PU_fastqc.html",path_in+"Ocotea-complicata-WE528_2PU_fastqc.html",path_in+"Ocotea-javitensis-WE529_2PU_fastqc.html",path_in+"Ocotea-skutchii-WE530_2PU_fastqc.html",path_in+"Ocotea-sinuata-WE531_2PU_fastqc.html",path_in+"Ocotea-botrantha-WE532_2PU_fastqc.html",path_in+"Nectandra-lineatifolia-WE533_2PU_fastqc.html",path_in+"Ocotea-foetens-WE521_UN_fastqc.html",path_in+"Ocotea-gabonensis-WE522_UN_fastqc.html",path_in+"Ocotea-meziana-WE523_UN_fastqc.html",path_in+"Pleurothyrium-cuneifolium-WE524_UN_fastqc.html",path_in+"Mespilodaphne-cymbarum-WE525_UN_fastqc.html",path_in+"Damburneya-gentlei-WE526_UN_fastqc.html",path_in+"Ocotea-glaucosericea-WE527_UN_fastqc.html",path_in+"Ocotea-complicata-WE528_UN_fastqc.html",path_in+"Ocotea-javitensis-WE529_UN_fastqc.html",path_in+"Ocotea-skutchii-WE530_UN_fastqc.html",path_in+"Ocotea-sinuata-WE531_UN_fastqc.html",path_in+"Ocotea-botrantha-WE532_UN_fastqc.html",path_in+"Nectandra-lineatifolia-WE533_UN_fastqc.html"] 
    outputs = [path_out+"multiqc_report.html", done]
    options = {'cores': 1, 'memory': "8g", 'walltime': "00:30:00", 'account':"cryptocarya"}


    spec = """

    source /home/laurakf/miniconda3/etc/profile.d/conda.sh

    conda activate multiqc

    multiqc -o {path_out} {path_in}
    
    touch {done}

    """.format(path_in = path_in, path_out = path_out, done = done)

    return (inputs, outputs, options, spec)



########################################################################################################################
################################################---- Hybpiper ----######################################################
########################################################################################################################
def hybpiper(name, p1, p2, un, path_out, path_in, done):
    """Hybpiper."""
    path_ins = [path_in+name+p1, path_in+name+p2, path_in+name+un] # The files which the job will look for before it runs
    outputs = [path_out+name"/", done] # The files which will have to be created in order for the job to be "completed"
    options = {'cores': 2, 'memory': "8g", 'walltime': "10:00:00", 'account':"cryptocarya"} #Slurm commands

    spec = """

    source /home/laurakf/miniconda3/etc/profile.d/conda.sh

    conda activate HybPiper

    cd /scratch/$SLURM_JOBID

    cd $TMPDIR
        
    hybpiper assemble --cpu 2 --targetfile_dna /home/laurakf/cryptocarya/TargetFile/mega353.fasta --readfiles {p1} {p2} --unpaired {un} --prefix {name} --bwa

    cp --recursive --update {name} /home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/

    touch {done}
    touch {out}{name}

    """.format(name=name, p1=path_in+name+p1, p2=path_in+name+p2, un=path_in+name+un, out=path_out, done=done)


    return (path_ins, outputs, options, spec)

# ########################################################################################################################
# #############################################---- Paralogs ----#########################################################
# ########################################################################################################################

# def paralogs(species,path_in, done, no_paralogs, in_done):
#     """Find Paralog genes and write them in the file called paralog.txt"""
#     path_ins = [path_in + species, in_done]
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
#     
    
#     touch {done}

#     """.format(sp = species, done = done, path_in = path_in, np = no_paralogs)
#     return (path_ins, outputs, options, spec)

# def no_paralogs(species, path_in, done, no_paralogs):
#     """Wrapper script to continue pipeline when Hybpiper finds no paralogs"""
#     path_ins = [path_in + species]
#     outputs = [done]
#     options = {'cores': 2, 'memory': "10g", 'walltime': "0:05:00", 'account':"Coryphoideae"}

#     spec = """

#     touch {done}
#     touch {np}

#     """.format(done=done, np=no_paralogs)
#     return(path_ins, outputs, options, spec)

# # ########################################################################################################################
# # #############################################---- Intronerate ----######################################################
# # ########################################################################################################################

# def intronerate(species, path_in, done):
#     """Intronerate the sequences from hybpiper."""
#     path_ins = [path_in + species, path_in+"done/Hybpiper/"+species]
#     outputs = [done]
#     options = {'cores': 4, 'memory': "20g", 'walltime': "16:00:00", 'account':"Coryphoideae"}

#     spec = """
#     source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
#     conda activate base

#     cd {path_in}

#     python3 /home/owrisberg/Coryphoideae/github_code/HybPiper/intronerate.py --prefix {sp} &>> intronerate_out.txt
        
    
#     touch {done}

#     """.format(sp = species, done = done, path_in = path_in)

#     return (path_ins, outputs, options, spec)

# # ########################################################################################################################
# # #############################################---- Coverage ----#########################################################
# # ########################################################################################################################
# def coverage(species, path_in, path_out, done,all_bam,all_sorted_bam, all_sorted_bam_bai, bam, cov,fasta,fasta_amb,fasta_ann,fasta_bwt,fasta_pac,fasta_sa,trimmed_fasta,up_bam,dir_in,dir_out):
#     """Calculating coverage of sequences."""
#     path_ins = [path_in+species, path_in+"done/Intronerate/"+species]
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

#     return (path_ins, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################

#Species removed from pipeline as they had no gene recovery [3050,3188,3272,3300,3316, 3364, 3392, 3394]
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


#### Running multiqc on trimmed data
gwf.target_from_template('multiqc_trimmed', multiqc_trimmed(path_in= "/home/laurakf/cryptocarya/Workflow/Test/04_FastQC/slidingwindow/",
                                                    path_out = "/home/laurakf/cryptocarya/Workflow/Test/05_MultiQC/slidingwindow/",
                                                    done = "/home/laurakf/cryptocarya/Workflow/Test/05_MultiQC/slidingwindow/done/multiqc_trimmed"))


sp = ["Ocotea-foetens-WE521","Ocotea-gabonensis-WE522","Ocotea-meziana-WE523","Pleurothyrium-cuneifolium-WE524","Mespilodaphne-cymbarum-WE525","Damburneya-gentlei-WE526","Ocotea-glaucosericea-WE527","Ocotea-complicata-WE528","Ocotea-javitensis-WE529","Ocotea-skutchii-WE530","Ocotea-sinuata-WE531","Ocotea-botrantha-WE532","Nectandra-lineatifolia-WE533"] 


for i in range(len(sp)):
    
    #### Running Hybpiper
    gwf.target_from_template('Hybpiper_'+str(i), hybpiper(name = sp[i],
                                                        p1 = "_1P.fastq",
                                                        p2 = "_2P.fastq",
                                                        un = "_UN.fastq",
                                                        path_out= "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/",
                                                        path_in = "/home/laurakf/cryptocarya/Workflow/Test/03_Trimmomatic/slidingwindow/",
                                                        done = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/HybPiper/"+sp[i]))
                                                                      

    # #### Paralogs
    
    # gwf.target_from_template('Paralogs_'+sp[i], paralogs(name = sp[i],
    #                                                     path_in = "/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
    #                                                     done = "/home/laurakf/cryptocarya/Workflow/Test/06_HybPiper/done/Paralogs/"+sp[i],
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

