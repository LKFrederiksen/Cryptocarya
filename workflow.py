
# ----------------------------------------------------------------------------------------------------------------------------
# This workflow is used to transform the raw sequence data into sequences ready for alignment.
# Workflow is the following

# 1: FastQC quality check of the sequences with slidingwindow or (maxinfo)
# 2: MultiQC wrapper to combine FastQC reports in one html file
# 3: Trimming of the sequences using trimmomatic
# 4: FastQC quality check of the trimmed sequences
# 5: MultiQC wrapper to combine FastQC reports in one html file
# 6: Hybpiper, using assemble (including intronerate) to get exons and introns for each sample. Retrieved stats. Paralog retriever to check for paralogs.
# 7: Trimming for Coverage of sequencing and joining the exon and intron of each species to create supercontigs.
# 8: Using script sample2genes (pastataster) to retrieve the gene files including all sample sequences found for each gene.
# 9: MAFFT to align sequences
# 10: Trimal to trim alignments
# 11: Optrimal to trim alignments
# 12: CIalign to trim alignments
# 13: Taper to trim alignments
# 14: IQtree to generate a hypothesis for each gene and generate gene trees
# 15: Generate ASTRAL tree from gene trees

#
# Workflow to reconstruct the outgroup of the PAFTOL tree using supercontigs. 
# ----------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------------
# This code is a variant of species_workflow.py by Oscar Balslev Wrisberg and workflow_completed.py by Paola de Lima Ferreira 
# Edited by Laura Kragh Frederiksen 12/10/2022
# ----------------------------------------------------------------------------------------------------------------------------

from os import O_SYNC, name
from gwf import Workflow
import os.path
import math
import glob

gwf = Workflow()

# ########################################################################################################################
# ###############################################---- Fastqc quality check raw ----#######################################
# ########################################################################################################################
# def fastqc_raw(name,path_in ,path_out, done,):
#     """Quality checking using fastqc as this should work on individual species"""
#     path_ins = [path_in+name+"_R1.fastq", path_in+name+"_R2.fastq"]
#     outputs = [path_out+name+"_R1_fastqc.html",path_out+name+"_R2_fastqc.html", done]
#     options = {'cores': 1, 'memory': "10g", 'walltime': "01:00:00", 'account':"cryptocarya"}


#     spec = """

#     echo {name}

#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh

#     conda activate fastqc

#     fastqc -o {path_out} {path_in}{name}_R1.fastq {path_in}{name}_R2.fastq
    
#     echo touching {done}
    
#     touch {done}

#     """.format(path_in = path_in, name = name, path_out = path_out, done = done)

#     return (path_ins, outputs, options, spec)

# #######################################################################################################################
# #############################################---- Multiqc quality check raw ----#######################################
# #######################################################################################################################
# def multiqc_raw(path_in ,path_out, done,):
#     """Quality checking using multiqc"""
#     inputs = [path_in+"Alse-petio-PAFTOL_R1_fastqc.html", path_in+"Athe-mosch-PAFTOL_R1_fastqc.html", path_in+"Beil-pendu-PAFTOL_R1_fastqc.html", path_in+"Beil-tsang-PAFTOL_R1_fastqc.html", path_in+"Cary-tonki-PAFTOL_R1_fastqc.html", path_in+"Caly-flori-PAFTOL_R1_fastqc.html", path_in+"Cass-filif-PAFTOL_R1_fastqc.html", path_in+"Cinn-camph-PAFTOL_R1_fastqc.html", path_in+"Cryp-alba-PAFTOL_R1_fastqc.html", path_in+"Deha-haina-PAFTOL_R1_fastqc.html", path_in+"Endi-macro-PAFTOL_R1_fastqc.html", path_in+"Gomo-keule-PAFTOL_R1_fastqc.html", path_in+"Hern-nymph-PAFTOL_R1_fastqc.html", path_in+"Idio-austr-PAFTOL_R1_fastqc.html", path_in+"Laur-nobil-PAFTOL_R1_fastqc.html", path_in+"Mach-salic-PAFTOL_R1_fastqc.html", path_in+"Magn-grand-PAFTOL_R1_fastqc.html", path_in+"Mezi-ita-uba-PAFTOL_R1_fastqc.html", path_in+"Moll-gilgi-PAFTOL_R1_fastqc.html", path_in+"Moni-rotun-PAFTOL_R1_fastqc.html", path_in+"Myri-fragr-PAFTOL_R1_fastqc.html", path_in+"Neoc-cauda-PAFTOL_R1_fastqc.html", path_in+"Noth-umbel-PAFTOL_R1_fastqc.html", path_in+"Pers-borbo-PAFTOL_R1_fastqc.html", path_in+"Peum-boldu-PAFTOL_R1_fastqc.html", path_in+"Phoe-lance-PAFTOL_R1_fastqc.html", path_in+"Sipa-guian-PAFTOL_R1_fastqc.html", path_in+"Spar-botoc-PAFTOL_R1_fastqc.html", path_in+"Synd-chine-PAFTOL_R1_fastqc.html", path_in+"Tamb-ficus-PAFTOL_R1_fastqc.html", path_in+"Alse-petio-PAFTOL_R2_fastqc.html", path_in+"Athe-mosch-PAFTOL_R2_fastqc.html", path_in+"Beil-pendu-PAFTOL_R2_fastqc.html", path_in+"Beil-tsang-PAFTOL_R2_fastqc.html", path_in+"Cary-tonki-PAFTOL_R2_fastqc.html", path_in+"Caly-flori-PAFTOL_R2_fastqc.html", path_in+"Cass-filif-PAFTOL_R2_fastqc.html", path_in+"Cinn-camph-PAFTOL_R2_fastqc.html", path_in+"Cryp-alba-PAFTOL_R2_fastqc.html", path_in+"Deha-haina-PAFTOL_R2_fastqc.html", path_in+"Endi-macro-PAFTOL_R2_fastqc.html", path_in+"Gomo-keule-PAFTOL_R2_fastqc.html", path_in+"Hern-nymph-PAFTOL_R2_fastqc.html", path_in+"Idio-austr-PAFTOL_R2_fastqc.html", path_in+"Laur-nobil-PAFTOL_R2_fastqc.html", path_in+"Mach-salic-PAFTOL_R2_fastqc.html", path_in+"Magn-grand-PAFTOL_R2_fastqc.html", path_in+"Mezi-ita-uba-PAFTOL_R2_fastqc.html", path_in+"Moll-gilgi-PAFTOL_R2_fastqc.html", path_in+"Moni-rotun-PAFTOL_R2_fastqc.html", path_in+"Myri-fragr-PAFTOL_R2_fastqc.html", path_in+"Neoc-cauda-PAFTOL_R2_fastqc.html", path_in+"Noth-umbel-PAFTOL_R2_fastqc.html", path_in+"Pers-borbo-PAFTOL_R2_fastqc.html", path_in+"Peum-boldu-PAFTOL_R2_fastqc.html", path_in+"Phoe-lance-PAFTOL_R2_fastqc.html", path_in+"Sipa-guian-PAFTOL_R2_fastqc.html", path_in+"Spar-botoc-PAFTOL_R2_fastqc.html", path_in+"Synd-chine-PAFTOL_R2_fastqc.html", path_in+"Tamb-ficus-PAFTOL_R2_fastqc.html"]
#     outputs = [path_out+"multiqc_report.html", path_out+"multiqc_data/", done]
#     options = {'cores': 1, 'memory': "10g", 'walltime': "01:00:00", 'account':"cryptocarya"}


#     spec = """

#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh

#     conda activate multiqc

#     multiqc -o {path_out} {path_in}
    
#     echo touching {done}

#     touch {done}

#     """.format(path_in = path_in, path_out = path_out, done = done)

#     return (inputs, outputs, options, spec)


# #######################################################################################################################################
# ################################################---- Trimmomatic - slidingwindow----###################################################
# #######################################################################################################################################
# def trimmomatic(name, path_in, path_out, done):
#     """Trimming raw data using trimmomatic with custom adapter.
#     Afterwards combines paired and unpaired reads for forward and reverse reads respectively for each species 
#     to enable post-trimming secapr quality_check for comparability before and after trimming """
#     path_ins = [path_in+name+"_R1.fastq", path_in+name+"_R2.fastq"]
#     outputs = [path_out+name+"_UN.fastq",path_out+"secapr_postrim/"+name+"_UN.fastq", done, path_out+name+"_1P.fastq", path_out+name+"_2P.fastq"]
#     options = {'cores': 8, 'memory': "8g", 'walltime': "02:00:00", 'account':"cryptocarya"}

#     spec = """
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate trimmomatic

#     trimmomatic PE -threads 16 -phred33 {path_in}_R1.fastq {path_in}_R2.fastq -baseout {output}.fastq\
#     ILLUMINACLIP:/home/laurakf/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:1:30:7:1:true\
#     LEADING:30\
#     SLIDINGWINDOW:4:30\
#     MINLEN:40\
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

#     echo touching {done}
    
#     touch {done}

#     """.format(path_in = path_in+name, output = path_out+name, done = done, name = name, path_out = path_out)

#     return (path_ins, outputs, options, spec)


# # #######################################################################################################################################
# # ################################################---- Trimmomatic - Maxinfo ----########################################################
# # #######################################################################################################################################
# # def trimmomatic_maxinfo(name, path_in, path_out, done):
# #     """Trimming raw data using trimmomatic with custom adapter.
# #     Afterwards combines paired and unpaired reads for forward and reverse reads respectively for each species 
# #     to enable post-trimming secapr quality_check for comparability before and after trimming """
# #     path_ins = [path_in+name+"_R1.fastq", path_in+name+"_1PU.fastq"]
# #     outputs = [path_out+name+"_UN.fastq",path_out+"secapr_postrim/"+name+"_UN.fastq", done, path_out+name+"_1P.fastq", path_out+name+"_2P.fastq"]
# #     options = {'cores': 8, 'memory': "5g", 'walltime': "00:30:00", 'account':"cryptocarya"}

# #     spec = """
# #     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
# #     conda activate trimmomatic

# #     trimmomatic PE -threads 16 -phred33 {path_in}_R1.fastq {path_in}_1PU.fastq -baseout {output}.fastq\
# #     ILLUMINACLIP:/home/laurakf/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:1:30:7:1:true\
# #     LEADING:3\
# #     TRAILING:3\
# #     MAXINFO:40:0.8\
# #     MINLEN:36\
# #     2>> stderr_trim_loop_output.txt

# #     echo combining {path_out}{name}_1P.fastq and {path_out}{name}_1U.fastq into {path_out}secapr_postrim/{name}_1PU.fastq 
# #     cat {path_out}{name}_1P.fastq {path_out}{name}_1U.fastq > {path_out}secapr_postrim/{name}_1PU.fastq 

# #     echo combining {path_out}{name}_2P.fastq and {path_out}{name}_2U.fastq into {path_out}secapr_postrim/{name}_2PU.fastq 
# #     cat {path_out}{name}_2P.fastq {path_out}{name}_2U.fastq > {path_out}secapr_postrim/{name}_2PU.fastq

# #     echo combining {path_out}{name}_1U.fastq {path_out}{name}_2U.fastq > {path_out}{name}_UN.fastq
# #     cat {path_out}{name}_1U.fastq {path_out}{name}_2U.fastq > {path_out}{name}_UN.fastq
# #     cp {path_out}{name}_UN.fastq {path_out}secapr_postrim/


# #     echo Removing {path_out}{name}_1U.fastq
# #     rm {path_out}{name}_1U.fastq

# #     echo Removing {path_out}{name}_2U.fastq
# #     rm {path_out}{name}_2U.fastq

# #     echo touching {done}

# #     touch {done}

# #     """.format(path_in = path_in+name, output = path_out+name, done = done, name = name, path_out = path_out)

# #     return (path_ins, outputs, options, spec)


# ########################################################################################################################################
# ###########################################---- Fastqc quality check trimmed (slidingwindow) ----#######################################
# ########################################################################################################################################
# def fastqc_trimmed(name,path_in ,path_out, done,):
#      """Quality checking using fastqc as this should work on individual species"""
#      path_ins = [path_in+name+"_UN.fastq", path_in+name+"_1PU.fastq", path_in+name+"_2PU.fastq"] # The files gwf looks for before it runs.
#      outputs = [path_out+name+"_1PU_fastqc.html", path_out+name+"_2PU_fastqc.html",path_out+name+"_UN_fastqc.html", done]
#      options = {'cores': 1, 'memory': "10g", 'walltime': "01:00:00", 'account':"cryptocarya"}


#      spec = """

#      echo {name}
     
#      source /home/laurakf/miniconda3/etc/profile.d/conda.sh

#      conda activate fastqc

#      fastqc -o {output} {path_in}{name}_1PU.fastq {path_in}{name}_2PU.fastq {path_in}{name}_UN.fastq
    
#      echo touching {done}

#      touch {done}

#      """.format(path_in = path_in, name = name, output = path_out, done = done)

#      return (path_ins, outputs, options, spec)


# ##################################################################################################################################
# ###########################################---- Fastqc quality check trimmed (maxinfo) ----#######################################
# ##################################################################################################################################
# # def fastqc_trimmed_maxinfo(name,path_in ,path_out, done,):
# #      """Quality checking using fastqc as this should work on individual species"""
# #      path_ins = [path_in+name+"_UN.fastq", path_in+name+"_1PU.fastq", path_in+name+"_2PU.fastq"] # The files gwf looks for before it runs.
# #      outputs = [path_out+name+"_1PU_fastqc.html", path_out+name+"_2PU_fastqc.html",path_out+name+"_UN_fastqc.html", done]
# #      options = {'cores': 1, 'memory': "10g", 'walltime': "00:30:00", 'account':"cryptocarya"}


# #      spec = """

# #      echo {name}
     
# #      source /home/laurakf/miniconda3/etc/profile.d/conda.sh

# #      conda activate fastqc

# #      fastqc -o {output} {path_in}{name}_1PU.fastq {path_in}{name}_2PU.fastq {path_in}{name}_UN.fastq
    
# #      echo touching {done}

# #      touch {done}

# #      """.format(path_in = path_in, name = name, output = path_out, done = done)

# #      return (path_ins, outputs, options, spec)

# #######################################################################################################################################
# #########################################---- Multiqc quality check trimmed (slidingwindow) ----#######################################
# #######################################################################################################################################
# def multiqc_trimmed(path_in ,path_out, done):
#     """Quality checking using multiqc"""
#     inputs = [path_in+"Alse-petio-PAFTOL_trimmed.fasta", path_in+"Athe-mosch-PAFTOL_trimmed.fasta", path_in+"Beil-pendu-PAFTOL_trimmed.fasta", path_in+"Beil-tsang-PAFTOL_trimmed.fasta", path_in+"Cary-tonki-PAFTOL_trimmed.fasta", path_in+"Caly-flori-PAFTOL_trimmed.fasta", path_in+"Cass-filif-PAFTOL_trimmed.fasta", path_in+"Cinn-camph-PAFTOL_trimmed.fasta", path_in+"Cryp-alba-PAFTOL_trimmed.fasta", path_in+"Deha-haina-PAFTOL_trimmed.fasta", path_in+"Endi-macro-PAFTOL_trimmed.fasta", path_in+"Gomo-keule-PAFTOL_trimmed.fasta", path_in+"Hern-nymph-PAFTOL_trimmed.fasta", path_in+"Idio-austr-PAFTOL_trimmed.fasta", path_in+"Laur-nobil-PAFTOL_trimmed.fasta", path_in+"Mach-salic-PAFTOL_trimmed.fasta", path_in+"Magn-grand-PAFTOL_trimmed.fasta", path_in+"Mezi-ita-uba-PAFTOL_trimmed.fasta", path_in+"Moll-gilgi-PAFTOL_trimmed.fasta", path_in+"Moni-rotun-PAFTOL_trimmed.fasta", path_in+"Myri-fragr-PAFTOL_trimmed.fasta", path_in+"Neoc-cauda-PAFTOL_trimmed.fasta", path_in+"Noth-umbel-PAFTOL_trimmed.fasta", path_in+"Pers-borbo-PAFTOL_trimmed.fasta", path_in+"Peum-boldu-PAFTOL_trimmed.fasta", path_in+"Phoe-lance-PAFTOL_trimmed.fasta", path_in+"Sipa-guian-PAFTOL_trimmed.fasta", path_in+"Spar-botoc-PAFTOL_trimmed.fasta", path_in+"Synd-chine-PAFTOL_trimmed.fasta", path_in+"Tamb-ficus-PAFTOL_trimmed.fasta", path_in+"Alse-petio-PAFTOL_1PU_fastqc.html", path_in+"Athe-mosch-PAFTOL_1PU_fastqc.html", path_in+"Beil-pendu-PAFTOL_1PU_fastqc.html", path_in+"Beil-tsang-PAFTOL_1PU_fastqc.html", path_in+"Cary-tonki-PAFTOL_1PU_fastqc.html", path_in+"Caly-flori-PAFTOL_1PU_fastqc.html", path_in+"Cass-filif-PAFTOL_1PU_fastqc.html", path_in+"Cinn-camph-PAFTOL_1PU_fastqc.html", path_in+"Cryp-alba-PAFTOL_1PU_fastqc.html", path_in+"Deha-haina-PAFTOL_1PU_fastqc.html", path_in+"Endi-macro-PAFTOL_1PU_fastqc.html", path_in+"Gomo-keule-PAFTOL_1PU_fastqc.html", path_in+"Hern-nymph-PAFTOL_1PU_fastqc.html", path_in+"Idio-austr-PAFTOL_1PU_fastqc.html", path_in+"Laur-nobil-PAFTOL_1PU_fastqc.html", path_in+"Mach-salic-PAFTOL_1PU_fastqc.html", path_in+"Magn-grand-PAFTOL_1PU_fastqc.html", path_in+"Mezi-ita-uba-PAFTOL_1PU_fastqc.html", path_in+"Moll-gilgi-PAFTOL_1PU_fastqc.html", path_in+"Moni-rotun-PAFTOL_1PU_fastqc.html", path_in+"Myri-fragr-PAFTOL_1PU_fastqc.html", path_in+"Neoc-cauda-PAFTOL_1PU_fastqc.html", path_in+"Noth-umbel-PAFTOL_1PU_fastqc.html", path_in+"Pers-borbo-PAFTOL_1PU_fastqc.html", path_in+"Peum-boldu-PAFTOL_1PU_fastqc.html", path_in+"Phoe-lance-PAFTOL_1PU_fastqc.html", path_in+"Sipa-guian-PAFTOL_1PU_fastqc.html", path_in+"Spar-botoc-PAFTOL_1PU_fastqc.html", path_in+"Synd-chine-PAFTOL_1PU_fastqc.html", path_in+"Tamb-ficus-PAFTOL_1PU_fastqc.html", path_in+"Alse-petio-PAFTOL_2PU_fastqc.html", path_in+"Athe-mosch-PAFTOL_2PU_fastqc.html", path_in+"Beil-pendu-PAFTOL_2PU_fastqc.html", path_in+"Beil-tsang-PAFTOL_2PU_fastqc.html", path_in+"Cary-tonki-PAFTOL_2PU_fastqc.html", path_in+"Caly-flori-PAFTOL_2PU_fastqc.html", path_in+"Cass-filif-PAFTOL_2PU_fastqc.html", path_in+"Cinn-camph-PAFTOL_2PU_fastqc.html", path_in+"Cryp-alba-PAFTOL_2PU_fastqc.html", path_in+"Deha-haina-PAFTOL_2PU_fastqc.html", path_in+"Endi-macro-PAFTOL_2PU_fastqc.html", path_in+"Gomo-keule-PAFTOL_2PU_fastqc.html", path_in+"Hern-nymph-PAFTOL_2PU_fastqc.html", path_in+"Idio-austr-PAFTOL_2PU_fastqc.html", path_in+"Laur-nobil-PAFTOL_2PU_fastqc.html", path_in+"Mach-salic-PAFTOL_2PU_fastqc.html", path_in+"Magn-grand-PAFTOL_2PU_fastqc.html", path_in+"Mezi-ita-uba-PAFTOL_2PU_fastqc.html", path_in+"Moll-gilgi-PAFTOL_2PU_fastqc.html", path_in+"Moni-rotun-PAFTOL_2PU_fastqc.html", path_in+"Myri-fragr-PAFTOL_2PU_fastqc.html", path_in+"Neoc-cauda-PAFTOL_2PU_fastqc.html", path_in+"Noth-umbel-PAFTOL_2PU_fastqc.html", path_in+"Pers-borbo-PAFTOL_2PU_fastqc.html", path_in+"Peum-boldu-PAFTOL_2PU_fastqc.html", path_in+"Phoe-lance-PAFTOL_2PU_fastqc.html", path_in+"Sipa-guian-PAFTOL_2PU_fastqc.html", path_in+"Spar-botoc-PAFTOL_2PU_fastqc.html", path_in+"Synd-chine-PAFTOL_2PU_fastqc.html", path_in+"Tamb-ficus-PAFTOL_2PU_fastqc.html"]
#     outputs = [path_out+"multiqc_report.html", path_out+"multiqc_data/", done]
#     options = {'cores': 1, 'memory': "10g", 'walltime': "01:00:00", 'account':"cryptocarya"}


#     spec = """

#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh

#     conda activate multiqc

#     multiqc -o {path_out} {path_in}

#     echo touching {done}

#     touch {done}

#     """.format(path_in = path_in, path_out = path_out, done = done)

#     return (inputs, outputs, options, spec)


# # ########################################################################################################################################
# # ##########################################---- Multiqc quality check trimmed (maxinfo) ----#############################################
# # ########################################################################################################################################
# # def multiqc_trimmed_maxinfo(path_in ,path_out, done,):
# #     """Quality checking using multiqc"""
# #     inputs = [path_in+"Aspi-fungi-686AL1_1PU_fastqc.html", path_in+"Aspi-parvi-687AL1_1PU_fastqc.html", path_in+"Beil-appen-688AL1_1PU_fastqc.html", path_in+"Beil-berte-689AL1_1PU_fastqc.html", path_in+"Beil-brach-690AL1_1PU_fastqc.html", path_in+"Beil-brene-691AL1_1PU_fastqc.html", path_in+"Beil-dicty-692AL1_1PU_fastqc.html", path_in+"Beil-emarg-693AL1_1PU_fastqc.html", path_in+"Beil-fasci-694AL1_1PU_fastqc.html", path_in+"Beil-fulva-695AL1_1PU_fastqc.html", path_in+"Beil-furfu-696AL1_1PU_fastqc.html", path_in+"Beil-hengh-697AL1_1PU_fastqc.html", path_in+"Beil-latif-699AL1_1PU_fastqc.html", path_in+"Beil-latif-700AL1_1PU_fastqc.html", path_in+"Beil-linha-701AL1_1PU_fastqc.html", path_in+"Beil-linoc-702AL1_1PU_fastqc.html", path_in+"Beil-macro-703AL1_1PU_fastqc.html", path_in+"Beil-madag-704AL1_1PU_fastqc.html", path_in+"Beil-manni-705AL1_1PU_fastqc.html", path_in+"Beil-manni-706AL1_1PU_fastqc.html", path_in+"Beil-manni-707AL1_1PU_fastqc.html", path_in+"Beil-miers-708AL1_1PU_fastqc.html", path_in+"Beil-morat-709AL1_1PU_fastqc.html", path_in+"Beil-pauci-710AL1_1PU_fastqc.html", path_in+"Beil-pedic-711AL1_1PU_fastqc.html", path_in+"Beil-per-C-713AL1_1PU_fastqc.html", path_in+"Beil-perco-712AL1_1PU_fastqc.html", path_in+"Beil-purpu-714AL1_1PU_fastqc.html", path_in+"Beil-robus-715AL1_1PU_fastqc.html", path_in+"Beil-roxbu-716AL1_1PU_fastqc.html", path_in+"Beil-roxbu-717AL1_1PU_fastqc.html", path_in+"Beil-rufoh-718AL1_1PU_fastqc.html", path_in+"Beil-rugos-719AL1_1PU_fastqc.html", path_in+"Beil-sary-720AL1_1PU_fastqc.html", path_in+"Beil-seric-721AL1_1PU_fastqc.html", path_in+"Beil-tarai-722AL1_1PU_fastqc.html", path_in+"Beil-tarai-723AL1_1PU_fastqc.html", path_in+"Beil-tawa-724AL1_1PU_fastqc.html", path_in+"Beil-tawa-742AL1_1PU_fastqc.html", path_in+"Beil-tawar-725AL1_1PU_fastqc.html", path_in+"Beil-tilar-726AL1_1PU_fastqc.html", path_in+"Beil-tungf-727AL1_1PU_fastqc.html", path_in+"Beil-ugand-728AL1_1PU_fastqc.html", path_in+"Beil-velut-729AL1_1PU_fastqc.html", path_in+"Beil-volck-730AL1_1PU_fastqc.html", path_in+"Beil-yunna-731AL1_1PU_fastqc.html", path_in+"Cryp-acuti-743AL1_1PU_fastqc.html", path_in+"Cryp-alba-800AL1_1PU_fastqc.html", path_in+"Cryp-albi-745AL1_1PU_fastqc.html", path_in+"Cryp-ampl-746AL1_1PU_fastqc.html", path_in+"Cryp-asche-747AL1_1PU_fastqc.html", path_in+"Cryp-botel-748AL1_1PU_fastqc.html", path_in+"Cryp-calci-749AL1_1PU_fastqc.html", path_in+"Cryp-chine-750AL1_1PU_fastqc.html", path_in+"Cryp-citri-751AL1_1PU_fastqc.html", path_in+"Cryp-conci-752AL1_1PU_fastqc.html", path_in+"Cryp-densi-753AL1_1PU_fastqc.html", path_in+"Cryp-ferre-754AL1_1PU_fastqc.html", path_in+"Cryp-fusca-755AL1_1PU_fastqc.html", path_in+"Cryp-haina-756AL1_1PU_fastqc.html", path_in+"Cryp-horne-757AL1_1PU_fastqc.html", path_in+"Cryp-krame-758AL1_1PU_fastqc.html", path_in+"Cryp-lepto-759AL1_1PU_fastqc.html", path_in+"Cryp-liebe-732AL1_1PU_fastqc.html", path_in+"Cryp-litor-733AL1_1PU_fastqc.html", path_in+"Cryp-litor-734AL1_1PU_fastqc.html", path_in+"Cryp-mandi-735AL1_1PU_fastqc.html", path_in+"Cryp-medic-736AL1_1PU_fastqc.html", path_in+"Cryp-micra-738AL1_1PU_fastqc.html", path_in+"Cryp-mosch-739AL1_1PU_fastqc.html", path_in+"Cryp-niten-740AL1_1PU_fastqc.html", path_in+"Cryp-oubat-741AL1_1PU_fastqc.html", path_in+"Cryp-ovali-778AL1_1PU_fastqc.html", path_in+"Cryp-pauci-779AL1_1PU_fastqc.html", path_in+"Cryp-pervi-780AL1_1PU_fastqc.html", path_in+"Cryp-pervi-781AL1_1PU_fastqc.html", path_in+"Cryp-polyn-782AL1_1PU_fastqc.html", path_in+"Cryp-polyn-783AL1_1PU_fastqc.html", path_in+"Cryp-rhodo-784AL1_1PU_fastqc.html", path_in+"Cryp-riede-785AL1_1PU_fastqc.html", path_in+"Cryp-rigid-786AL1_1PU_fastqc.html", path_in+"Cryp-rolle-787AL1_1PU_fastqc.html", path_in+"Cryp-salig-788AL1_1PU_fastqc.html", path_in+"Cryp-sello-789AL1_1PU_fastqc.html", path_in+"Cryp-spath-790AL1_1PU_fastqc.html", path_in+"Cryp-spath-791AL1_1PU_fastqc.html", path_in+"Cryp-subtr-793AL1_1PU_fastqc.html", path_in+"Cryp-thou-794AL1_1PU_fastqc.html", path_in+"Cryp-trans-795AL1_1PU_fastqc.html", path_in+"Cryp-vello-796AL1_1PU_fastqc.html", path_in+"Cryp-woodi-797AL1_1PU_fastqc.html", path_in+"Cryp-wylie-798AL1_1PU_fastqc.html", path_in+"Cryp-yunna-799AL1_1PU_fastqc.html", path_in+"Endi-impre-801AL1_1PU_fastqc.html", path_in+"Endi-jones-802AL1_1PU_fastqc.html", path_in+"Endi-latif-803AL1_1PU_fastqc.html", path_in+"Endi-lecar-804AL1_1PU_fastqc.html", path_in+"Endi-palme-805AL1_1PU_fastqc.html", path_in+"Endi-phaeo-767AL1_1PU_fastqc.html", path_in+"Endi-pilos-760AL1_1PU_fastqc.html", path_in+"Endi-poueb-761AL1_1PU_fastqc.html", path_in+"Endi-puben-762AL1_1PU_fastqc.html", path_in+"Endi-sanke-763AL1_1PU_fastqc.html", path_in+"Endi-scrob-764AL1_1PU_fastqc.html", path_in+"Endi-sulav-765AL1_1PU_fastqc.html", path_in+"Endi-xanth-766AL1_1PU_fastqc.html", path_in+"Eusi-zwage-768AL1_1PU_fastqc.html", path_in+"Pota-confl-769AL1_1PU_fastqc.html", path_in+"Pota-micro-770AL1_1PU_fastqc.html", path_in+"Pota-obtus-771AL1_1PU_fastqc.html", path_in+"Pota-obtus-772AL1_1PU_fastqc.html", path_in+"Poto-melag-773AL1_1PU_fastqc.html", path_in+"Sino-hongk-774AL1_1PU_fastqc.html", path_in+"Synd-kwang-775AL1_1PU_fastqc.html", path_in+"Synd-marli-776AL1_1PU_fastqc.html", path_in+"Synd-marli-777AL1_1PU_fastqc.html", path_in+"Aspi-fungi-686AL1_2PU_fastqc.html", path_in+"Aspi-parvi-687AL1_2PU_fastqc.html", path_in+"Beil-appen-688AL1_2PU_fastqc.html", path_in+"Beil-berte-689AL1_2PU_fastqc.html", path_in+"Beil-brach-690AL1_2PU_fastqc.html", path_in+"Beil-brene-691AL1_2PU_fastqc.html", path_in+"Beil-dicty-692AL1_2PU_fastqc.html", path_in+"Beil-emarg-693AL1_2PU_fastqc.html", path_in+"Beil-fasci-694AL1_2PU_fastqc.html", path_in+"Beil-fulva-695AL1_2PU_fastqc.html", path_in+"Beil-furfu-696AL1_2PU_fastqc.html", path_in+"Beil-hengh-697AL1_2PU_fastqc.html", path_in+"Beil-latif-699AL1_2PU_fastqc.html", path_in+"Beil-latif-700AL1_2PU_fastqc.html", path_in+"Beil-linha-701AL1_2PU_fastqc.html", path_in+"Beil-linoc-702AL1_2PU_fastqc.html", path_in+"Beil-macro-703AL1_2PU_fastqc.html", path_in+"Beil-madag-704AL1_2PU_fastqc.html", path_in+"Beil-manni-705AL1_2PU_fastqc.html", path_in+"Beil-manni-706AL1_2PU_fastqc.html", path_in+"Beil-manni-707AL1_2PU_fastqc.html", path_in+"Beil-miers-708AL1_2PU_fastqc.html", path_in+"Beil-morat-709AL1_2PU_fastqc.html", path_in+"Beil-pauci-710AL1_2PU_fastqc.html", path_in+"Beil-pedic-711AL1_2PU_fastqc.html", path_in+"Beil-per-C-713AL1_2PU_fastqc.html", path_in+"Beil-perco-712AL1_2PU_fastqc.html", path_in+"Beil-purpu-714AL1_2PU_fastqc.html", path_in+"Beil-robus-715AL1_2PU_fastqc.html", path_in+"Beil-roxbu-716AL1_2PU_fastqc.html", path_in+"Beil-roxbu-717AL1_2PU_fastqc.html", path_in+"Beil-rufoh-718AL1_2PU_fastqc.html", path_in+"Beil-rugos-719AL1_2PU_fastqc.html", path_in+"Beil-sary-720AL1_2PU_fastqc.html", path_in+"Beil-seric-721AL1_2PU_fastqc.html", path_in+"Beil-tarai-722AL1_2PU_fastqc.html", path_in+"Beil-tarai-723AL1_2PU_fastqc.html", path_in+"Beil-tawa-724AL1_2PU_fastqc.html", path_in+"Beil-tawa-742AL1_2PU_fastqc.html", path_in+"Beil-tawar-725AL1_2PU_fastqc.html", path_in+"Beil-tilar-726AL1_2PU_fastqc.html", path_in+"Beil-tungf-727AL1_2PU_fastqc.html", path_in+"Beil-ugand-728AL1_2PU_fastqc.html", path_in+"Beil-velut-729AL1_2PU_fastqc.html", path_in+"Beil-volck-730AL1_2PU_fastqc.html", path_in+"Beil-yunna-731AL1_2PU_fastqc.html", path_in+"Cryp-acuti-743AL1_2PU_fastqc.html", path_in+"Cryp-alba-800AL1_2PU_fastqc.html", path_in+"Cryp-albi-745AL1_2PU_fastqc.html", path_in+"Cryp-ampl-746AL1_2PU_fastqc.html", path_in+"Cryp-asche-747AL1_2PU_fastqc.html", path_in+"Cryp-botel-748AL1_2PU_fastqc.html", path_in+"Cryp-calci-749AL1_2PU_fastqc.html", path_in+"Cryp-chine-750AL1_2PU_fastqc.html", path_in+"Cryp-citri-751AL1_2PU_fastqc.html", path_in+"Cryp-conci-752AL1_2PU_fastqc.html", path_in+"Cryp-densi-753AL1_2PU_fastqc.html", path_in+"Cryp-ferre-754AL1_2PU_fastqc.html", path_in+"Cryp-fusca-755AL1_2PU_fastqc.html", path_in+"Cryp-haina-756AL1_2PU_fastqc.html", path_in+"Cryp-horne-757AL1_2PU_fastqc.html", path_in+"Cryp-krame-758AL1_2PU_fastqc.html", path_in+"Cryp-lepto-759AL1_2PU_fastqc.html", path_in+"Cryp-liebe-732AL1_2PU_fastqc.html", path_in+"Cryp-litor-733AL1_2PU_fastqc.html", path_in+"Cryp-litor-734AL1_2PU_fastqc.html", path_in+"Cryp-mandi-735AL1_2PU_fastqc.html", path_in+"Cryp-medic-736AL1_2PU_fastqc.html", path_in+"Cryp-micra-738AL1_2PU_fastqc.html", path_in+"Cryp-mosch-739AL1_2PU_fastqc.html", path_in+"Cryp-niten-740AL1_2PU_fastqc.html", path_in+"Cryp-oubat-741AL1_2PU_fastqc.html", path_in+"Cryp-ovali-778AL1_2PU_fastqc.html", path_in+"Cryp-pauci-779AL1_2PU_fastqc.html", path_in+"Cryp-pervi-780AL1_2PU_fastqc.html", path_in+"Cryp-pervi-781AL1_2PU_fastqc.html", path_in+"Cryp-polyn-782AL1_2PU_fastqc.html", path_in+"Cryp-polyn-783AL1_2PU_fastqc.html", path_in+"Cryp-rhodo-784AL1_2PU_fastqc.html", path_in+"Cryp-riede-785AL1_2PU_fastqc.html", path_in+"Cryp-rigid-786AL1_2PU_fastqc.html", path_in+"Cryp-rolle-787AL1_2PU_fastqc.html", path_in+"Cryp-salig-788AL1_2PU_fastqc.html", path_in+"Cryp-sello-789AL1_2PU_fastqc.html", path_in+"Cryp-spath-790AL1_2PU_fastqc.html", path_in+"Cryp-spath-791AL1_2PU_fastqc.html", path_in+"Cryp-subtr-793AL1_2PU_fastqc.html", path_in+"Cryp-thou-794AL1_2PU_fastqc.html", path_in+"Cryp-trans-795AL1_2PU_fastqc.html", path_in+"Cryp-vello-796AL1_2PU_fastqc.html", path_in+"Cryp-woodi-797AL1_2PU_fastqc.html", path_in+"Cryp-wylie-798AL1_2PU_fastqc.html", path_in+"Cryp-yunna-799AL1_2PU_fastqc.html", path_in+"Endi-impre-801AL1_2PU_fastqc.html", path_in+"Endi-jones-802AL1_2PU_fastqc.html", path_in+"Endi-latif-803AL1_2PU_fastqc.html", path_in+"Endi-lecar-804AL1_2PU_fastqc.html", path_in+"Endi-palme-805AL1_2PU_fastqc.html", path_in+"Endi-phaeo-767AL1_2PU_fastqc.html", path_in+"Endi-pilos-760AL1_2PU_fastqc.html", path_in+"Endi-poueb-761AL1_2PU_fastqc.html", path_in+"Endi-puben-762AL1_2PU_fastqc.html", path_in+"Endi-sanke-763AL1_2PU_fastqc.html", path_in+"Endi-scrob-764AL1_2PU_fastqc.html", path_in+"Endi-sulav-765AL1_2PU_fastqc.html", path_in+"Endi-xanth-766AL1_2PU_fastqc.html", path_in+"Eusi-zwage-768AL1_2PU_fastqc.html", path_in+"Pota-confl-769AL1_2PU_fastqc.html", path_in+"Pota-micro-770AL1_2PU_fastqc.html", path_in+"Pota-obtus-771AL1_2PU_fastqc.html", path_in+"Pota-obtus-772AL1_2PU_fastqc.html", path_in+"Poto-melag-773AL1_2PU_fastqc.html", path_in+"Sino-hongk-774AL1_2PU_fastqc.html", path_in+"Synd-kwang-775AL1_2PU_fastqc.html", path_in+"Synd-marli-776AL1_2PU_fastqc.html", path_in+"Synd-marli-777AL1_2PU_fastqc.html", path_in+"Aspi-fungi-686AL1_UN_fastqc.html", path_in+"Aspi-parvi-687AL1_UN_fastqc.html", path_in+"Beil-appen-688AL1_UN_fastqc.html", path_in+"Beil-berte-689AL1_UN_fastqc.html", path_in+"Beil-brach-690AL1_UN_fastqc.html", path_in+"Beil-brene-691AL1_UN_fastqc.html", path_in+"Beil-dicty-692AL1_UN_fastqc.html", path_in+"Beil-emarg-693AL1_UN_fastqc.html", path_in+"Beil-fasci-694AL1_UN_fastqc.html", path_in+"Beil-fulva-695AL1_UN_fastqc.html", path_in+"Beil-furfu-696AL1_UN_fastqc.html", path_in+"Beil-hengh-697AL1_UN_fastqc.html", path_in+"Beil-latif-699AL1_UN_fastqc.html", path_in+"Beil-latif-700AL1_UN_fastqc.html", path_in+"Beil-linha-701AL1_UN_fastqc.html", path_in+"Beil-linoc-702AL1_UN_fastqc.html", path_in+"Beil-macro-703AL1_UN_fastqc.html", path_in+"Beil-madag-704AL1_UN_fastqc.html", path_in+"Beil-manni-705AL1_UN_fastqc.html", path_in+"Beil-manni-706AL1_UN_fastqc.html", path_in+"Beil-manni-707AL1_UN_fastqc.html", path_in+"Beil-miers-708AL1_UN_fastqc.html", path_in+"Beil-morat-709AL1_UN_fastqc.html", path_in+"Beil-pauci-710AL1_UN_fastqc.html", path_in+"Beil-pedic-711AL1_UN_fastqc.html", path_in+"Beil-per-C-713AL1_UN_fastqc.html", path_in+"Beil-perco-712AL1_UN_fastqc.html", path_in+"Beil-purpu-714AL1_UN_fastqc.html", path_in+"Beil-robus-715AL1_UN_fastqc.html", path_in+"Beil-roxbu-716AL1_UN_fastqc.html", path_in+"Beil-roxbu-717AL1_UN_fastqc.html", path_in+"Beil-rufoh-718AL1_UN_fastqc.html", path_in+"Beil-rugos-719AL1_UN_fastqc.html", path_in+"Beil-sary-720AL1_UN_fastqc.html", path_in+"Beil-seric-721AL1_UN_fastqc.html", path_in+"Beil-tarai-722AL1_UN_fastqc.html", path_in+"Beil-tarai-723AL1_UN_fastqc.html", path_in+"Beil-tawa-724AL1_UN_fastqc.html", path_in+"Beil-tawa-742AL1_UN_fastqc.html", path_in+"Beil-tawar-725AL1_UN_fastqc.html", path_in+"Beil-tilar-726AL1_UN_fastqc.html", path_in+"Beil-tungf-727AL1_UN_fastqc.html", path_in+"Beil-ugand-728AL1_UN_fastqc.html", path_in+"Beil-velut-729AL1_UN_fastqc.html", path_in+"Beil-volck-730AL1_UN_fastqc.html", path_in+"Beil-yunna-731AL1_UN_fastqc.html", path_in+"Cryp-acuti-743AL1_UN_fastqc.html", path_in+"Cryp-alba-800AL1_UN_fastqc.html", path_in+"Cryp-albi-745AL1_UN_fastqc.html", path_in+"Cryp-ampl-746AL1_UN_fastqc.html", path_in+"Cryp-asche-747AL1_UN_fastqc.html", path_in+"Cryp-botel-748AL1_UN_fastqc.html", path_in+"Cryp-calci-749AL1_UN_fastqc.html", path_in+"Cryp-chine-750AL1_UN_fastqc.html", path_in+"Cryp-citri-751AL1_UN_fastqc.html", path_in+"Cryp-conci-752AL1_UN_fastqc.html", path_in+"Cryp-densi-753AL1_UN_fastqc.html", path_in+"Cryp-ferre-754AL1_UN_fastqc.html", path_in+"Cryp-fusca-755AL1_UN_fastqc.html", path_in+"Cryp-haina-756AL1_UN_fastqc.html", path_in+"Cryp-horne-757AL1_UN_fastqc.html", path_in+"Cryp-krame-758AL1_UN_fastqc.html", path_in+"Cryp-lepto-759AL1_UN_fastqc.html", path_in+"Cryp-liebe-732AL1_UN_fastqc.html", path_in+"Cryp-litor-733AL1_UN_fastqc.html", path_in+"Cryp-litor-734AL1_UN_fastqc.html", path_in+"Cryp-mandi-735AL1_UN_fastqc.html", path_in+"Cryp-medic-736AL1_UN_fastqc.html", path_in+"Cryp-micra-738AL1_UN_fastqc.html", path_in+"Cryp-mosch-739AL1_UN_fastqc.html", path_in+"Cryp-niten-740AL1_UN_fastqc.html", path_in+"Cryp-oubat-741AL1_UN_fastqc.html", path_in+"Cryp-ovali-778AL1_UN_fastqc.html", path_in+"Cryp-pauci-779AL1_UN_fastqc.html", path_in+"Cryp-pervi-780AL1_UN_fastqc.html", path_in+"Cryp-pervi-781AL1_UN_fastqc.html", path_in+"Cryp-polyn-782AL1_UN_fastqc.html", path_in+"Cryp-polyn-783AL1_UN_fastqc.html", path_in+"Cryp-rhodo-784AL1_UN_fastqc.html", path_in+"Cryp-riede-785AL1_UN_fastqc.html", path_in+"Cryp-rigid-786AL1_UN_fastqc.html", path_in+"Cryp-rolle-787AL1_UN_fastqc.html", path_in+"Cryp-salig-788AL1_UN_fastqc.html", path_in+"Cryp-sello-789AL1_UN_fastqc.html", path_in+"Cryp-spath-790AL1_UN_fastqc.html", path_in+"Cryp-spath-791AL1_UN_fastqc.html", path_in+"Cryp-subtr-793AL1_UN_fastqc.html", path_in+"Cryp-thou-794AL1_UN_fastqc.html", path_in+"Cryp-trans-795AL1_UN_fastqc.html", path_in+"Cryp-vello-796AL1_UN_fastqc.html", path_in+"Cryp-woodi-797AL1_UN_fastqc.html", path_in+"Cryp-wylie-798AL1_UN_fastqc.html", path_in+"Cryp-yunna-799AL1_UN_fastqc.html", path_in+"Endi-impre-801AL1_UN_fastqc.html", path_in+"Endi-jones-802AL1_UN_fastqc.html", path_in+"Endi-latif-803AL1_UN_fastqc.html", path_in+"Endi-lecar-804AL1_UN_fastqc.html", path_in+"Endi-palme-805AL1_UN_fastqc.html", path_in+"Endi-phaeo-767AL1_UN_fastqc.html", path_in+"Endi-pilos-760AL1_UN_fastqc.html", path_in+"Endi-poueb-761AL1_UN_fastqc.html", path_in+"Endi-puben-762AL1_UN_fastqc.html", path_in+"Endi-sanke-763AL1_UN_fastqc.html", path_in+"Endi-scrob-764AL1_UN_fastqc.html", path_in+"Endi-sulav-765AL1_UN_fastqc.html", path_in+"Endi-xanth-766AL1_UN_fastqc.html", path_in+"Eusi-zwage-768AL1_UN_fastqc.html", path_in+"Pota-confl-769AL1_UN_fastqc.html", path_in+"Pota-micro-770AL1_UN_fastqc.html", path_in+"Pota-obtus-771AL1_UN_fastqc.html", path_in+"Pota-obtus-772AL1_UN_fastqc.html", path_in+"Poto-melag-773AL1_UN_fastqc.html", path_in+"Sino-hongk-774AL1_UN_fastqc.html", path_in+"Synd-kwang-775AL1_UN_fastqc.html", path_in+"Synd-marli-776AL1_UN_fastqc.html", path_in+"Synd-marli-777AL1_UN_fastqc.html"] 
# #     outputs = [path_out+"multiqc_report.html", path_out+"multiqc_data/", done]
# #     options = {'cores': 1, 'memory': "10g", 'walltime': "00:30:00", 'account':"cryptocarya"}


# #     spec = """

# #     source /home/laurakf/miniconda3/etc/profile.d/conda.sh

# #     conda activate multiqc

# #     multiqc -o {path_out} {path_in}

# #     echo touching {done}

# #     touch {done}

# #     """.format(path_in = path_in, path_out = path_out, done = done)

# #     return (inputs, outputs, options, spec)



# ########################################################################################################################
# ################################################---- Hybpiper ----######################################################
# ########################################################################################################################
# def hybpiper(name, p1, p2, un, path_out, path_in, done):
#     """Hybpiper."""
#     path_ins = [path_in+name+p1, path_in+name+p2, path_in+name+un] # The files which the job will look for before it runs
#     outputs = [path_out+name, done] # The files which will have to be created in order for the job to be "completed"
#     options = {'cores': 2, 'memory': "12g", 'walltime': "2:00:00", 'account':"cryptocarya"} #Slurm commands

#     spec = """

#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh

#     conda activate HybPiper

#     TMPDIR=/scratch/$SLURM_JOBID
#     export TMPDIR
#     mkdir -p $TMPDIR
#     cd $TMPDIR
    
#     # Here I have used the Rohwer target file!
#     hybpiper assemble --cpu 2 --targetfile_dna /home/laurakf/cryptocarya/TargetFile/mega353_rohwer.fasta --readfiles {p1} {p2} --unpaired {un} --prefix {name} --bwa --run_intronerate

#     cp --recursive --update {name} /home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/

#     echo touching {done}

#     touch {done}
    

#     """.format(name=name, p1=path_in+name+p1, p2=path_in+name+p2, un=path_in+name+un, out=path_out+name, done=done)


#     return (path_ins, outputs, options, spec)

# ########################################################################################################################
# ###################################################---- Stats ----######################################################
# ########################################################################################################################

# # In this step you should run the statistics on the folder where we have the Hybpiper_results
# # I did not create a folder just for Hybpiper results, then I will create here and move the assemble results to there

# def stats(path_in, done, path_out, in_done, name):
#    """Gather statistics about the HybPiper run(s).""", 
#    path_ins = [path_in+name, in_done] # The files that has to be present before the job runs.
#    outputs = [path_out+"seq_lengths.tsv", path_out+"hybpiper_stats.tsv", path_out+"recovery_heatmap.png"]  # The files which will have to be created in order for the job to be "completed"
#    options = {'cores': 2, 'memory': "16g", 'walltime': "04:00:00", 'account':"cryptocarya"} #Slurm commands

#    spec = """
   
#    source /home/laurakf/miniconda3/etc/profile.d/conda.sh

#    conda activate HybPiper
    
#    cd {path_in}
    
#    hybpiper stats --targetfile_dna /home/laurakf/cryptocarya/TargetFile/mega353_rohwer.fasta 'supercontig' {path_in}namelist.txt # Get stats

#    hybpiper recovery_heatmap {path_in}seq_lengths.tsv # Make heatmap

#    mv seq_lengths.tsv {path_out} # Move all stats and the heatmap to a new subfolder
    
#    mv hybpiper_stats.tsv {path_out}

#    mv recovery_heatmap.png {path_out} 

#    echo touching {done}

#    touch {done}
      
#    """.format(path_in = path_in, done = done, path_out = path_out, in_done = in_done, name = name)

#    return (path_ins, outputs, options, spec) 


# ########################################################################################################################
# #############################################---- Paralogs ----#########################################################
# ########################################################################################################################

# ##### This is the approach to use. #####

# def paralogs(name, path_in, done, in_done, path_out):
#     """Run HybPiper v. 2.1 - paralog retriever """
#     path_ins = [path_in+name, in_done]
#     outputs = [done]
#     options = {'cores': 2, 'memory': "10g", 'walltime': "01:00:00", 'account':"cryptocarya"}

#     spec = """
    
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
    
#     conda activate HybPiper
    
#     cd {path_in}

#     hybpiper paralog_retriever namelist.txt -t_dna /home/laurakf/cryptocarya/TargetFile/mega353_rohwer.fasta
    
#     mv paralog_report.tsv {path_out}
#     mv paralogs_above_threshold_report.txt {path_out}
#     mv paralogs_all {path_out}
#     mv paralogs_no_chimeras {path_out}
#     mv paralog_heatmap.png {path_out}

#     echo touching {done}

#     touch {done}

#      """.format(name = name, done = done, path_in = path_in, path_out = path_out, in_done = in_done)
    
#     return (path_ins, outputs, options, spec)

# ### Genes found with paralogs: 4951, 4989, 5343, 5347, 5355, 5428, 5434, 5859, 5940, 5958, 6110, 6387, 6449, 6498, 6782, 6955, 6995, 7336.


# ########################################################################################################################
# #############################################---- Coverage ----#########################################################
# ########################################################################################################################

# #This script does the following:
# # Gather all contigs from each sample in one fasta file: coverage/sample.fasta
# # Map paired and unpaired reads to that fasta using BWA mem
# # Deduplicate reads using Picard
# # Calculate depth using samtools
# # Mask/strip any bases with coverage <2
# # Generate a new trimmed sample-level fasta: coverage/sample_trimmed.fasta

# def coverage(name, path_in, path_out, done,all_bam,all_sorted_bam, all_sorted_bam_bai, bam, cov,fasta,fasta_amb,fasta_ann,fasta_bwt,fasta_pac,fasta_sa,trimmed_fasta,up_bam,dir_in,dir_out, dir_wrk):
#     """Calculating coverage of sequences."""
#     path_ins = [path_in+name]
#     outputs = [path_out+name+all_bam,
#      path_out+name+all_sorted_bam,
#       path_out+name+all_sorted_bam_bai,
#        path_out+name+bam,
#     path_out+name+cov,
#      path_out+name+fasta,
#       path_out+name+fasta_amb,
#        path_out+name+fasta_ann,
#         path_out+name+fasta_bwt,
#     path_out+name+fasta_pac,
#      path_out+name+fasta_sa,
#       path_out+name+trimmed_fasta,
#        path_out+name+up_bam,done] #ALL the output files
#     options = {'cores': 4, 'memory': "24g", 'walltime': "01:30:00", 'account':"cryptocarya"}

#     spec = """
    
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
    
#     conda activate HybPiper
    
#     cd {path_in}

#     python3 /home/laurakf/cryptocarya/Scripts/coverage.py {name} {dir_in} {dir_out} {dir_wrk}
    
#     echo touching {done}

#     touch {done}

#     """.format(name = name, done = done, path_in = path_in, dir_in = dir_in, dir_out = dir_out, dir_wrk = dir_wrk)

#     return (path_ins, outputs, options, spec)


# ########################################################################################################################
# #############################################---- Retrieve ----#########################################################
# ########################################################################################################################

# #Think about doing blacklisting here? you could just remove species from the inputs here if you dont want them in the downstream analysis
# # I will add Synd-chine, Peum-boldu and Chim-salic as gene files from PAFTOL before retrieving. Look at Chim-salic to see an example.

# def retrieve(path_in, done):
#     """Retrieve gene sequences from all the species and create an unaligned multifasta for each gene."""
#     path_ins = [path_in+"Alse-petio-PAFTOL_trimmed.fasta", path_in+"Athe-mosch-PAFTOL_trimmed.fasta", path_in+"Beil-pendu-PAFTOL_trimmed.fasta", path_in+"Beil-tsang-PAFTOL_trimmed.fasta", path_in+"Cary-tonki-PAFTOL_trimmed.fasta", path_in+"Caly-flori-PAFTOL_trimmed.fasta", path_in+"Cass-filif-PAFTOL_trimmed.fasta", path_in+"Chim-salic-PAFTOL_trimmed.fasta", path_in+"Cinn-camph-PAFTOL_trimmed.fasta", path_in+"Cryp-alba-PAFTOL_trimmed.fasta", path_in+"Deha-haina-PAFTOL_trimmed.fasta", path_in+"Endi-macro-PAFTOL_trimmed.fasta", path_in+"Gomo-keule-PAFTOL_trimmed.fasta", path_in+"Hern-nymph-PAFTOL_trimmed.fasta", path_in+"Idio-austr-PAFTOL_trimmed.fasta", path_in+"Laur-nobil-PAFTOL_trimmed.fasta", path_in+"Mach-salic-PAFTOL_trimmed.fasta", path_in+"Magn-grand-PAFTOL_trimmed.fasta", path_in+"Mezi-ita-uba-PAFTOL_trimmed.fasta", path_in+"Moll-gilgi-PAFTOL_trimmed.fasta", path_in+"Moni-rotun-PAFTOL_trimmed.fasta", path_in+"Myri-fragr-PAFTOL_trimmed.fasta", path_in+"Neoc-cauda-PAFTOL_trimmed.fasta", path_in+"Noth-umbel-PAFTOL_trimmed.fasta", path_in+"Pers-borbo-PAFTOL_trimmed.fasta", path_in+"Peum-boldu-PAFTOL_trimmed.fasta", path_in+"Phoe-lance-PAFTOL_trimmed.fasta", path_in+"Sipa-guian-PAFTOL_trimmed.fasta", path_in+"Spar-botoc-PAFTOL_trimmed.fasta", path_in+"Synd-chine-PAFTOL_trimmed.fasta", path_in+"Tamb-ficus-PAFTOL_trimmed.fasta"]
#     outputs = [done]
#     options = {'cores': 4, 'memory': "8g", 'walltime': "00:30:00", 'account':"cryptocarya"}

#     spec = """
    
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh

#     conda activate HybPiper

#     cd {path_in}

#     ls *trimmed.fasta > filelist.txt

#     python3 /home/laurakf/cryptocarya/Scripts/sample2genes.py > outstats.csv

#     echo touching {done}

#     touch {done}

#     """.format(path_in = path_in, done = done)

#     return (path_ins, outputs, options, spec)
    
# I have included Chim-salic, Synd-chine and Peum-boldu from PAFTOL gene sequences. (No UN sequence which causes an error when running HybPiper).
### Here you should wait for the output. The output will comprise a file for each gene with the species sequence recovered.


# ##########################################################################################################################
# ###############################################---- MAFFT ----#############################################################
# ##########################################################################################################################

# # Here go to folder 08_Retrieve and ls -1
# # Get the gene names and write them in genes = []
# # We found  genes for Lauraceae using the mega353 target file.

# def mafft(genes, path_in, path_out, done, gene):
#     """Aligning all the sequences for each gene."""
#     path_ins = [path_in+genes]
#     outputs = [done, path_out+gene+"_aligned.fasta"] 
#     options = {'cores': 4, 'memory': "4g", 'walltime': "1:00:00", 'account':"cryptocarya"}

#     spec = """

#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh

#     conda activate Mafft

#     cd {path_in}

#     mafft --thread 4 --globalpair --adjustdirectionaccurately --maxiterate 1000 {genes} > {path_out}{gene}_aligned.fasta

#     echo touching {done}

#     touch {done}

#     """.format(genes = genes, done = done, path_in = path_in, path_out = path_out, gene = gene)

#     return (path_ins, outputs, options, spec)
    
# #It is a good idea to rename the fasta files here.


# ########################################################################################################################
# ###############################################---- TRIMAL ----#########################################################
# ########################################################################################################################

# #Cleaning according trimal
# #Get raw alignments and trim them according to a gap threshold.

# # Make sure to adjust path in script 'gap_trimming.sh'
# # Seems like it is necessary to rename the *.fasta to *.fasta.old yourself. Otherwise perhaps my path was wrong?
# # find . -depth -name "*.fasta", -exec sh -c 'f="{}"; mv -- "$f", "${f%.fasta}.fasta.old"' \;
# # This way it did take the *.old files and move them to each directory 0.1 etc. and trim them. 

# def gt_trimming(path_in, path_out, done, gene):
#     """, Use trimal for trimming all alignments for each of the GT values specified"""
#     inputs = [path_in+gene+"_aligned.fasta"]
#     outputs = [done]
#     options = {'cores': 1, 'memory': "5g", 'walltime': "0:20:00", 'account':"cryptocarya"}

#     spec="""
    
#     #Activating Trimal
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate Trimal

#     # Moved alignments to Trimal folder
#     #Go to alignments folder 
#     cd {path_in}

#     #Running gaptrimming.sh
#     cd {path_out}
#     bash /home/laurakf/cryptocarya/Scripts/gap_trimming.sh -g {gene}_aligned.fasta.old

#     echo touching {done}

#     touch {done}

#     """.format(path_in = path_in, done = done, gene = gene, path_out = path_out)

#     return(inputs, outputs, options, spec)
    
    
#######################################################################################################################
################################################---- AMAS ----#########################################################
#######################################################################################################################
#Make file summary_0.txt

# ####### Calculating amas summary (raw)

# #For raw alignments
# def amas_raw(path_in, done, in_done):
#     """Creating summary files for all the trimmed alignments for each raw alignment"""
#     inputs = [in_done]
#     outputs = [path_in+"summary_0.txt", done]
#     options = {'cores': 1, 'memory': "2g", 'walltime': "0:10:00", 'account':"cryptocarya"}

#     spec="""

#     #Activating AMAS
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate Amas
    
#     cd {path_in}

#     #Calculating amas summary
#     /home/laurakf/cryptocarya/Scripts/AMAS/amas/AMAS.py summary -f fasta -d dna -i *.fasta.old

#     mv summary.txt summary_0.txt 

#     echo touching {done}

#     touch {done}
    
#     """.format(path_in = path_in, done = done, in_done = in_done)

#     return(inputs, outputs, options, spec)


# Removed   5599_aligned.fasta.old,  5921_aligned.fasta.old, 5943_aligned.fasta.old, 6041_aligned.fasta.old, 
#Hereafter you need to do some manual work and remove the headlines statistics.


# ######## Calculating amas summary (gt)
# # Make file summary_0.1.txt, summary_0.15.txt, summary_0.2.txt etc. 

# #For cutoff (gt) alignments
# def amas_gt(path_in, cut_off, done, in_done):
#     """Creating summary files for all the trimmed alignments for each raw alignment"""
#     inputs = [path_in+cut_off, in_done]
#     outputs = [path_in+"summary_"+cut_off+".txt", done]
#     options = {'cores': 1, 'memory': "2g", 'walltime': "0:10:00", 'account':"cryptocarya"}

#     spec="""

#     #Activating AMAS
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate Amas
    
#     cd {path_in}{cut_off} 

#     #Calculating amas summary
#     /home/laurakf/cryptocarya/Scripts/AMAS/amas/AMAS.py summary -f fasta -d dna -i *.fasta.old
   
#     mv summary.txt ../summary_{cut_off}.txt 

#     touch {done}

#     echo {done}

#     echo {cut_off}

#     echo {path_in}
    
#     """.format(path_in = path_in, cut_off = cut_off, done = done, in_done = in_done)

#     return(inputs, outputs, options, spec)


# ########################################################################################################################
# #############################################---- Optrimal ----#########################################################
# ########################################################################################################################

# # Run the script optrimal.R in interactive mode from the folder 10_trimal where the summary files generated in the 2 amas scripts are.
# # It is important that the cutoff.txt file has the labels 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9
# # and that the summary files are named summary_0.txt, summary_0.1.txt, summary_0.15.txt, summary_0.2.txt etc. 
# # Make folder called optimal_final_results.
# # The interactive mode is called by activating the conda environment R: conda activate R. When write R --interactive. Hereafter insert all the lines from the optrimal.R script.
# # Make done file in folder optrimal.


# #Getting the best alignment for each gene 
# def optrim(path_in, path_out, done):
#     """Select the best alignments according to the gt value"""
#     inputs = [path_in+"summary_0.txt", path_in+"summary_0.1.txt", path_in+"summary_0.15.txt",path_in+"summary_0.2.txt",path_in+"summary_0.25.txt",path_in+"summary_0.3.txt",
#     path_in+"summary_0.35.txt",path_in+"summary_0.4.txt",path_in+"summary_0.45.txt",path_in+"summary_0.5.txt",path_in+"summary_0.55.txt",path_in+"summary_0.6.txt",path_in+"summary_0.65.txt",
#     path_in+"summary_0.7.txt",path_in+"summary_0.75.txt",path_in+"summary_0.8.txt",path_in+"summary_0.85.txt",path_in+"summary_0.9.txt"]
#     outputs = [done, path_out+"optimal_final_results"]
#     options = {'cores': 2, 'memory': "5g", 'walltime': "0:10:00", 'account':"cryptocarya"}

#     spec="""

#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate base

#  #Going to folder with trimmed files
#     cd {path_in}

#     Rscript --vanilla /home/laurakf/cryptocarya/Scripts/optrimal.R

#     echo touching {done}

#     touch {done}

#     """.format(path_in = path_in, path_out = path_out, done = done)

#     return(inputs, outputs, options, spec)

#Here we should include remove empty files?

#HERE you can have an annoying error in R: "Error in real_loss[1:(length(real_loss) - 1)] :  only 0's may be mixed with negative subscripts"
# Download all the summary and cutoff files.
# Run the R code in your computer
# Check the file lost, and see if there is a gene that only has 1.00000 for all gt values, and delete those
# Copy the genes that you delete from raw_alignments to the trimal_output


# ####################################################################################################################################################
# ##############################################---- Move optrimal files to new folder ----###########################################################
# ####################################################################################################################################################

# def move(path_in, path_out, in_done, done):
#     """Moving files from trimal folder to optrimal folder."""
#     inputs = [in_done]
#     outputs = [done]
#     options = {'cores': 1, 'memory': "2g", 'walltime': "00:05:00", 'account':"cryptocarya"}

#     spec = """

#     cd {path_in}

#     mv dldp_* {path_out}optrim_output/
    
#     mv optimal_final_results {path_out}

#     mv overlost.txt {path_out}

#     echo touching {done}

#     touch {done}

#     """.format(path_in = path_in, path_out = path_out, done = done, in_done = in_done)

#     return(inputs, outputs, options, spec)

# ##########################################################################################################################
# ##############################################---- CIALIGN ----###########################################################
# ##########################################################################################################################

# def cialign1(gene, path_in, path_out, done):
#     """Cleaning alignments using cialign default."""
#     inputs = [path_in+gene+"_aligned.fasta.old"]
#     outputs = [path_out+gene+"_cialign_cleaned_cleaned.fasta", done]
#     options = {'cores': 4, 'memory': "5g", 'walltime': "00:15:00", 'account':"cryptocarya"}

#     spec = """
    
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate CIAlign

#     cd {path_in}

#     CIAlign --infile {gene}_aligned.fasta.old --all --outfile_stem {path_out}{gene}_cialign_cleaned

#     echo touching {done}

#     touch {done}

#     """.format(gene = gene, done = done, path_in = path_in, path_out = path_out)

#     return (inputs, outputs, options, spec)
    
# # Here we lost several genes.

# ########################################################################################################################
# ###############################################---- TAPER ----##########################################################
# ########################################################################################################################

# def taper(path_in, gene, path_out, done):
#     """Using TAPER AFTER CIAlign to remove errors in small species-specific stretches of the multiple sequence alignments"""
#     inputs = [path_in+gene+"_cialign_cleaned_cleaned.fasta"]
#     outputs = [done]
#     options = {'cores': 1, 'memory': "5g", 'walltime': "00:30:00", 'account':"cryptocarya"}

#     spec = """
     
#     cd {path_in}
    
#     # Activate Taper    
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate Taper
        
#     julia /home/laurakf/cryptocarya/Programs/TAPER-master/correction_multi.jl {gene}_cialign_cleaned_cleaned.fasta > {gene}_output_taper.fasta 
    
#     mv {gene}_output_taper.fasta {path_out}

#     echo touching {done}

#     touch {done}
        
#     """.format(path_in = path_in, gene = gene, path_out = path_out, done = done)

#     return (inputs, outputs, options, spec)
    

# ########################################################################################################################
# #############################################---- Exon Mapper ----######################################################
# ########################################################################################################################

# def exon_map(path_in,path_out,done,gene):
#     # """This creates new alignments in `07_mapping` that contain the original alignments plus the exon sequences of the two species that had the highest recovery success at each locus.."""
#     inputs = ["/home/laurakf/cryptocarya/Workflow/PAFTOL/13_Taper/done/"+gene]
#     outputs = [done,path_out+gene+"_output_taper_mapped.fasta"] 
#     options = {'cores': 4, 'memory': "10g", 'walltime': "01:00:00", 'account':"cryptocarya"}

#     spec="""

#     #Activating conda HybPiper environment (to get numpy etc.)
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate HybPiper

#     #Going to folder with data (Taper)
#     cd {path_in}

#     # Running Exon_mapper
#     python3 /home/laurakf/cryptocarya/Scripts/exon_mapper.py --gene {gene} --outdir {path_out} --file_ending _output_taper.fasta

#     touch {done}
#     """.format(path_in=path_in, done=done, gene=gene, path_out=path_out)

#     return(inputs, outputs, options, spec)


# ##########################################################################################################################
# ##########################---- Copying alignments and creating partitions ----############################################
# ##########################################################################################################################

# def partitioner(path_in, path_out, gene, done):
#     """Copying alignments from the manual alignment folder to the treebuilding folder and creating partition files"""
#     inputs = [path_in+gene+"_output_taper.fasta"]
#     outputs = [path_out+gene+"_part.txt",path_out+gene+"_clean.fasta", done]
#     options = {'cores': 1, 'memory': "5g", 'walltime': "00:20:00", 'account':"cryptocarya"}

#     spec = """

#         source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#         conda activate HybPiper


#         #Going to folder with data
#         cd {path_in}
    
#         #The partitioner should produce 2 files for each gene
#         #one file called gene_aligned_part.txt which is the partitioning file
#         #another called gene_aligned_clean.fasta which are just the sequences without the exons

#         python3 /home/laurakf/cryptocarya/Scripts/partitioner.py --smoother 10 --gene {gene} --file_ending _output_taper.fasta
        
#         #Moving files to the correct folders
#         mv {gene}_clean.fasta {path_out}
#         mv {gene}_part.txt {path_out} 

#         touch {done}


#     """.format(path_in = path_in, gene = gene, done = done, path_out = path_out)

#     return (inputs, outputs, options, spec)


# ###########################################################################################################################
# #############################################---- IQ-tree ----#############################################################
# ###########################################################################################################################

# def iq_tree(path_in, gene,path_out, done):
#     """Using Iq-tree to produce trees for each gene with a partition file to use individual substitution rates for each gene"""
#     inputs = [path_in+gene+"_part.txt", path_in+gene+"_clean.fasta"]
#     outputs = [path_out+gene+".txt.tre"]
#     options = {'cores': 10, 'memory': "10g", 'walltime': "01:00:00", 'account':"cryptocarya"}

#     spec = """


#     # Activate IQtree    
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate IQtree

#     cd {path_in}

#     #Actual IQtree tree search. 
#     iqtree2 -s {gene}_clean.fasta -p {gene}_part.txt -T AUTO -ntmax 20 -m MFP -B 1000 -redo 

#     mv {gene}*.treefile {path_out}
#     mv {gene}*.model.gz {path_out}
#     mv {gene}*.contree {path_out}
#     mv {gene}*.bionj {path_out}
#     mv {gene}*.ckp.gz {path_out}
#     mv {gene}*.iqtree {path_out}
#     mv {gene}*.log {path_out}
#     mv {gene}*.mldist {path_out}
#     mv {gene}*.splits.nex {path_out}
        
#     touch {done}

#     """.format(path_in = path_in, gene = gene, path_out=path_out, done = done)

#     return (inputs, outputs, options, spec)


# ########################################################################################################################
# #####################################---- Astral Tree Search ----#######################################################
# ########################################################################################################################


# def astral(path_in, done):
#     """Using Astral to construct a species tree based on the gene trees"""
#     inputs = [path_in]
#     outputs = [done]
#     options = {'cores': 20, 'memory': "40g", 'walltime': "1:00:00", 'account':"cryptocarya"}

#     spec = """
                     
#     cd {path_in} 
  
    # while read name
    # do cat RAxML_bipartitions.*_tree >> beilschmiedia_trees.tre && echo "", >> beilschmiedia_trees.tre
#     done < /home/laurakf/cryptocarya/Workflow/Lauraceae/15_Raxml/Astral_names.txt

#     while read name
#     do cat RAxML_bestTree.*_tree >> beilschmiedia_bestTree_trees.tre && echo "", >> beilschmiedia_bestTree_trees.tre
#     done < /home/laurakf/cryptocarya/Workflow/Lauraceae/15_Raxml/Astral_names.txt

#     /home/laurakf/cryptocarya/Programs/newick-utils-1.6/src/nw_ed beilschmiedia_trees.tre 'i & b<=10' o > beilschmiedia_trees_BP10.tre

#     java -jar /home/laurakf/cryptocarya/Programs/Astral/astral.5.7.8.jar -i beilschmiedia_trees_BP10.tre -t 2 -o beilschmiedia_trees_BP10_SpeciesTree_annotQ.tre

#     java -jar /home/laurakf/cryptocarya/Programs/Astral/astral.5.7.8.jar -i beilschmiedia_trees_BP10.tre -t 0 -o beilschmiedia_trees_BP10_SpeciesTree.tre

#     /home/laurakf/cryptocarya/Programs/phyx/src/pxrr -t beilschmiedia_trees_BP10_SpeciesTree.tre -g Myri-fragr-PAFTOL, Magn-grand-PAFTOL > beilschmiedia_trees_BP10_SpeciesTree_rooted.tre

#     /home/laurakf/cryptocarya/Programs/phyx/src/pxrr -t beilschmiedia_trees_BP10_SpeciesTree_annotQ.tre -g Myri-fragr-PAFTOL, Magn-grand-PAFTOL > beilschmiedia_trees_BP10_SpeciesTree_annotQ_rooted.tre


#     touch {done}

#     """.format(path_in = path_in, done = done)

#     return (inputs, outputs, options, spec)



##### Code beneath should be run to retrieve sortadate output. Can among others be used for tree calibration. #####

########################################################################################################################
#####################################---- Rooting gene trees----########################################################
########################################################################################################################
# def root_genetrees(path_in, gene, path_out):
#     """Using rerooter.py to root each individual gene tree based on the available outgroup"""
#     inputs = [path_in+gene+"_part.txt.treefile"] # changed _part.txt.tre to .txt.tre
#     outputs = [path_in+gene+"_rooted.tre"]
#     options = {'cores': 2, 'memory': "5g", 'walltime': "00:30:00", 'account':"cryptocarya"}

#     spec = """

#     # Activating Conda environment
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate Phyx


#     cd {path_in}

#     echo Rerooting each genetree based on the outgroup  
#     python3 /home/laurakf/cryptocarya/Scripts/rooter.py --gene {gene} --treefile {gene}_part.txt.treefile 

#     mv {gene}_rooted.tre {path_out}


#     """.format(path_in = path_in, gene = gene, path_out = path_out)

#     return (inputs, outputs, options, spec)


# #############################################################################################################################
# #####################################---- SortaDate ----#####################################################################
# #############################################################################################################################

# def sorta_date(path_in, path_out, astral_tree, done):
#     """Using SortaDate to produce a CSV file which can be used to evaluate the use of different genes in dating the trees"""
#     inputs = [astral_tree]
#     outputs = [path_out+"var",path_out+"bp",path_out+"comb", path_out+"gg", done]
#     options = {'cores': 3, 'memory': "10g", 'walltime': "00:10:00", 'account':"cryptocarya"}

#     spec = """

#     #Activating conda base environment 
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate Phyx # SortaDate scripts are dependent on python and 3 phyx programs. 

#     #Get the root-to-tip variance with
#     python {SortaDate}get_var_length.py {path_in} --flend _rooted.tre --outf {path_out}var --outg Magn-grand-PAFTOL,Myri-fragr-PAFTOL

#     #Get the bipartition support with
#     python {SortaDate}get_bp_genetrees.py {path_in} {astral_tree} --flend _rooted.tre --outf {path_out}bp

#     #Combine the results from these two runs with
#     python {SortaDate}combine_results.py {path_out}var {path_out}bp --outf {path_out}comb

#     #Sort and get the list of the good genes with
#     python {SortaDate}get_good_genes.py {path_out}comb --max 1000 --order 3,1,2 --outf {path_out}gg

#     touch {done}


#     """.format(path_in = path_in, path_out = path_out, SortaDate = "/home/laurakf/cryptocarya/Programs/SortaDate/src/", astral_tree = astral_tree, done = done)

#     return (inputs, outputs, options, spec)


# #######################--- Code to reconstruct outgroup with combination of supercontig partitioner and exon data ---#######################

# Need to run gene trees for exons in IQtree format to make the format of supercontig partitioner and exon gene trees match.

# #################################################################################################################################
# #############################################---- IQ-tree (exon)----#############################################################
# #################################################################################################################################

# def iq_tree_exon(path_in, gene,path_out, done):
#     """Using Iq-tree to produce trees for each gene with a partition file to use individual substitution rates for each gene"""
#     inputs = [path_in+gene+"_output_taper.fasta"]
#     outputs = [path_out+gene+"_output_taper.fasta.treefile"]
#     options = {'cores': 8, 'memory': "4g", 'walltime': "01:00:00", 'account':"cryptocarya"}

#     spec = """


#     # Activate IQtree    
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate IQtree

#     cd {path_in}

#     #Actual IQtree tree search. 
#     iqtree2 -s {gene}_output_taper.fasta -T AUTO -ntmax 20 -m MFP -B 1000 -redo 

#     mv {gene}*.treefile {path_out}
#     mv {gene}*.model.gz {path_out}
#     mv {gene}*.contree {path_out}
#     mv {gene}*.bionj {path_out}
#     mv {gene}*.ckp.gz {path_out}
#     mv {gene}*.iqtree {path_out}
#     mv {gene}*.log {path_out}
#     mv {gene}*.mldist {path_out}
#     mv {gene}*.splits.nex {path_out}
        
#     touch {done}

#     """.format(path_in = path_in, gene = gene, path_out=path_out, done = done)

#     return (inputs, outputs, options, spec)


# Prior to running the below code a folder with the required genes - a combination of supercontig partitioner and exon genes should be made.
# First try with a folder excluding the 10% worst supercontig partitioner genes and replace them with the exon genes. 
# Repeat until the correct topology has been obtained.

########################################################################################################################
##############################################---- Folder ----##########################################################
########################################################################################################################

# Folder: -10% supercontig partitioner genes
# Supercontig partitioner genes: 6139 6500 4951 5614 6036 7029 7111 5116 6064 4992 5162 6098 6488 5427 5335 6913 5131 5620 6318 5770 6363 5594 5168 5554 5328 6016 6563 6460 6384 7313 5858 6782 6420 5032 5318 6883 5333 7572 6494 5639 5339 4954 5899 6909 4802 5716 6373 5404 6462 6570 5271 5699 6528 5849 6639 6992 5463 5426 5893 5428 6914 4527 6961 6544 5449 5034 6051 7136 5702 7577 5454 6533 6457 5163 5894 6924 4932 5348 6447 7371 5421 5280 5870 6496 5791 5945 5634 7363 6072 5038 5343 5913 5596 7324 6458 5857 5816 5942 5981 6295 6641 4942 5815 5406 6689 5326 4757 6779 4793 6954 7024 6652 6848 6933 6003 6978 7279 5926 5502 6068 6404 5477 6738 6631 5138 6649 5977 6432 6450 6459 6216 6968 7067 7021 6376 5177 7141 7336 6865 7174 6732 4691 4724 5910 6947 5562 4848 6175 6282 6038 6412 6483 6667 5464 5264 5664 6860 6882 5865 5104 5670 5220 6733 6636 5703 5551 5950 6389 6559 5200 7241 5644 6572 5018 5721 7331 6454 5366 6685 6439 6825 5257 5528 7367 5296 5206 4890 6000 6298 6492 5918 5843 5273 5802 7325 7333 6299 6378 5536 5513 4889 6527 5933 6004 5968 6550 6284 6383 5531 6148 6226 5299 6274 6538 5744 5489 6176 6797 6029 5469 5188 5398 6601 6227 5949 5578 5980 5990 7194 6398 6050 5841 5840 6532 6198 6258 6303 5821 6552 6265 5960 6238 6320 6875 5919 6717 5090 6792 4796 7628 6979 6746 6407 5944 6854 5460 4806 6366 6405 6507 6048 5822 6946 
# Exons: 5866 6859 5123 6026 5974 6962 5304 6164 5772 6056 6785 5853 6114 6128 7128 6034 4471 6130 6660 6958 4744 7135 6506 6620 6221 6713 6393 5354 5660 7273 6526

# Folder: -20% supercontig partitioner genes
# Supercontig partitioner genes: 6139 6500 4951 5614 6036 7029 7111 5116 6064 4992 5162 6098 6488 5427 5335 6913 5131 5620 6318 5770 6363 5594 5168 5554 5328 6016 6563 6460 6384 7313 5858 6782 6420 5032 5318 6883 5333 7572 6494 5639 5339 4954 5899 6909 4802 5716 6373 5404 6462 6570 5271 5699 6528 5849 6639 6992 5463 5426 5893 5428 6914 4527 6961 6544 5449 5034 6051 7136 5702 7577 5454 6533 6457 5163 5894 6924 4932 5348 6447 7371 5421 5280 5870 6496 5791 5945 5634 7363 6072 5038 5343 5913 5596 7324 6458 5857 5816 5942 5981 6295 6641 4942 5815 5406 6689 5326 4757 6779 4793 6954 7024 6652 6848 6933 6003 6978 7279 5926 5502 6068 6404 5477 6738 6631 5138 6649 5977 6432 6450 6459 6216 6968 7067 7021 6376 5177 7141 7336 6865 7174 6732 4691 4724 5910 6947 5562 4848 6175 6282 6038 6412 6483 6667 5464 5264 5664 6860 6882 5865 5104 5670 5220 6733 6636 5703 5551 5950 6389 6559 5200 7241 5644 6572 5018 5721 7331 6454 5366 6685 6439 6825 5257 5528 7367 5296 5206 4890 6000 6298 6492 5918 5843 5273 5802 7325 7333 6299 6378 5536 5513 4889 6527 5933 6004 5968 6550 6284 6383 5531 6148 6226 5299 6274 6538 5744 5489 6176 6797 6029 5469 5188 5398 6601 6227 5949 5578 5980 5990 7194 6398 6050 5841 5840 6532 6198 6258 6303 5821 6552 6265 5960 6238 6320 6875 5919 6717 5090 6792 4796 7628 6979 6746 6407 5944 6854 5460 4806 6366 6405 6507 6048 5822 6946 
# Exons: 6198 6258 6303 5821 6552 6265 5960 6238 6320 6875 5919 6717 5090 6792 4796 7628 6979 6746 6407 5944 6854 5460 4806 6366 6405 6507 6048 5822 6946 5866 6859 5123 6026 5974 6962 5304 6164 5772 6056 6785 5853 6114 6128 7128 6034 4471 6130 6660 6958 4744 7135 6506 6620 6221 6713 6393 5354 5660 7273 6526

# Folder: -30% supercontig partitioner genes
# Supercontig partitioner genes: 6139 6500 4951 5614 6036 7029 7111 5116 6064 4992 5162 6098 6488 5427 5335 6913 5131 5620 6318 5770 6363 5594 5168 5554 5328 6016 6563 6460 6384 7313 5858 6782 6420 5032 5318 6883 5333 7572 6494 5639 5339 4954 5899 6909 4802 5716 6373 5404 6462 6570 5271 5699 6528 5849 6639 6992 5463 5426 5893 5428 6914 4527 6961 6544 5449 5034 6051 7136 5702 7577 5454 6533 6457 5163 5894 6924 4932 5348 6447 7371 5421 5280 5870 6496 5791 5945 5634 7363 6072 5038 5343 5913 5596 7324 6458 5857 5816 5942 5981 6295 6641 4942 5815 5406 6689 5326 4757 6779 4793 6954 7024 6652 6848 6933 6003 6978 7279 5926 5502 6068 6404 5477 6738 6631 5138 6649 5977 6432 6450 6459 6216 6968 7067 7021 6376 5177 7141 7336 6865 7174 6732 4691 4724 5910 6947 5562 4848 6175 6282 6038 6412 6483 6667 5464 5264 5664 6860 6882 5865 5104 5670 5220 6733 6636 5703 5551 5950 6389 6559 5200 7241 5644 6572 5018 5721 7331 6454 5366 6685 6439 6825 5257 5528 7367 5296 5206 4890 6000 6298 6492 5918 5843 5273 5802 7325 7333 6299 6378 5536 5513 4889 6527 5933 6004 5968 6550 6284 6383 5531 6148 6226 5299 6274 6538 5744 5489 6176 6797 6029 5469 5188 5398 6601 6227 5949 5578 5980 5990 7194 6398 6050 5841 5840 6532 6198 6258 6303 5821 6552 6265 5960 6238 6320 6875 5919 6717 5090 6792 4796 7628 6979 6746 6407 5944 6854 5460 4806 6366 6405 6507 6048 5822 6946 
# Exons: "5968", "6550", "6284", "6383", "5531", "6148", "6226", "5299", "6274", "6538", "5744", "5489", "6176", "6797", "6029", "5469", "5188", "5398", "6601", "6227", "5949", "5578", "5980", "5990", "7194", "6398", "6050", "5841", "5840", "6532", "6198", "6258", "6303", "5821", "6552", "6265", "5960", "6238", "6320", "6875", "5919", "6717", "5090", "6792", "4796", "7628", "6979", "6746", "6407", "5944", "6854", "5460", "4806", "6366", "6405", "6507", "6048", "5822", "6946", "5866", "6859", "5123", "6026", "5974", "6962", "5304", "6164", "5772", "6056", "6785", "5853", "6114", "6128", "7128", "6034", "4471", "6130", "6660", "6958", "4744", "7135", "6506", "6620", "6221", "6713", "6393", "5354", "5660", "7273", "6526"

# Folder: -40% supercontig partitioner genes
# Supercontig partitioner genes: 6139 6500 4951 5614 6036 7029 7111 5116 6064 4992 5162 6098 6488 5427 5335 6913 5131 5620 6318 5770 6363 5594 5168 5554 5328 6016 6563 6460 6384 7313 5858 6782 6420 5032 5318 6883 5333 7572 6494 5639 5339 4954 5899 6909 4802 5716 6373 5404 6462 6570 5271 5699 6528 5849 6639 6992 5463 5426 5893 5428 6914 4527 6961 6544 5449 5034 6051 7136 5702 7577 5454 6533 6457 5163 5894 6924 4932 5348 6447 7371 5421 5280 5870 6496 5791 5945 5634 7363 6072 5038 5343 5913 5596 7324 6458 5857 5816 5942 5981 6295 6641 4942 5815 5406 6689 5326 4757 6779 4793 6954 7024 6652 6848 6933 6003 6978 7279 5926 5502 6068 6404 5477 6738 6631 5138 6649 5977 6432 6450 6459 6216 6968 7067 7021 6376 5177 7141 7336 6865 7174 6732 4691 4724 5910 6947 5562 4848 6175 6282 6038 6412 6483 6667 5464 5264 5664 6860 6882 5865 5104 5670 5220 6733 6636 5703 5551 5950 6389 6559 5200 7241 5644 6572 5018 5721 7331 6454 5366 6685 6439 6825 5257 5528 7367 5296 5206 4890 6000 6298 6492 5918 5843 5273 5802 7325 7333 6299 6378 5536 5513 4889 6527 5933 6004 5968 6550 6284 6383 5531 6148 6226 5299 6274 6538 5744 5489 6176 6797 6029 5469 5188 5398 6601 6227 5949 5578 5980 5990 7194 6398 6050 5841 5840 6532 6198 6258 6303 5821 6552 6265 5960 6238 6320 6875 5919 6717 5090 6792 4796 7628 6979 6746 6407 5944 6854 5460 4806 6366 6405 6507 6048 5822 6946 
# Exons: "7331", "6454", "5366", "6685", "6439", "6825", "5257", "5528", "7367", "5296", "5206", "4890", "6000", "6298", "6492", "5918", "5843", "5273", "5802", "7325", "7333", "6299", "6378", "5536", "5513", "4889", "6527", "5933", "6004", "5968", "6550", "6284", "6383", "5531", "6148", "6226", "5299", "6274", "6538", "5744", "5489", "6176", "6797", "6029", "5469", "5188", "5398", "6601", "6227", "5949", "5578", "5980", "5990", "7194", "6398", "6050", "5841", "5840", "6532", "6198", "6258", "6303", "5821", "6552", "6265", "5960", "6238", "6320", "6875", "5919", "6717", "5090", "6792", "4796", "7628", "6979", "6746", "6407", "5944", "6854", "5460", "4806", "6366", "6405", "6507", "6048", "5822", "6946", "5866", "6859", "5123", "6026", "5974", "6962", "5304", "6164", "5772", "6056", "6785", "5853", "6114", "6128", "7128", "6034", "4471", "6130", "6660", "6958", "4744", "7135", "6506", "6620", "6221", "6713", "6393", "5354", "5660", "7273", "6526

# Folder: -50% supercontig partitioner genes
# Supercontig partitioner genes: 6139 6500 4951 5614 6036 7029 7111 5116 6064 4992 5162 6098 6488 5427 5335 6913 5131 5620 6318 5770 6363 5594 5168 5554 5328 6016 6563 6460 6384 7313 5858 6782 6420 5032 5318 6883 5333 7572 6494 5639 5339 4954 5899 6909 4802 5716 6373 5404 6462 6570 5271 5699 6528 5849 6639 6992 5463 5426 5893 5428 6914 4527 6961 6544 5449 5034 6051 7136 5702 7577 5454 6533 6457 5163 5894 6924 4932 5348 6447 7371 5421 5280 5870 6496 5791 5945 5634 7363 6072 5038 5343 5913 5596 7324 6458 5857 5816 5942 5981 6295 6641 4942 5815 5406 6689 5326 4757 6779 4793 6954 7024 6652 6848 6933 6003 6978 7279 5926 5502 6068 6404 5477 6738 6631 5138 6649 5977 6432 6450 6459 6216 6968 7067 7021 6376 5177 7141 7336 6865 7174 6732 4691 4724 5910 6947 5562 4848 6175 6282 6038 6412 6483 6667 5464 5264 5664 6860 6882 5865 5104 5670 5220 6733 6636 5703 5551 5950 6389 6559 5200 7241 5644 6572 5018 5721 7331 6454 5366 6685 6439 6825 5257 5528 7367 5296 5206 4890 6000 6298 6492 5918 5843 5273 5802 7325 7333 6299 6378 5536 5513 4889 6527 5933 6004 5968 6550 6284 6383 5531 6148 6226 5299 6274 6538 5744 5489 6176 6797 6029 5469 5188 5398 6601 6227 5949 5578 5980 5990 7194 6398 6050 5841 5840 6532 6198 6258 6303 5821 6552 6265 5960 6238 6320 6875 5919 6717 5090 6792 4796 7628 6979 6746 6407 5944 6854 5460 4806 6366 6405 6507 6048 5822 6946 
# Exons: "4848", "6175", "6282", "6038", "6412", "6483", "6667", "5464", "5264", "5664", "6860", "6882", "5865", "5104", "5670", "5220", "6733", "6636", "5703", "5551", "5950", "6389", "6559", "5200", "7241", "5644", "6572", "5018", "5721", "7331", "6454", "5366", "6685", "6439", "6825", "5257", "5528", "7367", "5296", "5206", "4890", "6000", "6298", "6492", "5918", "5843", "5273", "5802", "7325", "7333", "6299", "6378", "5536", "5513", "4889", "6527", "5933", "6004", "5968", "6550", "6284", "6383", "5531", "6148", "6226", "5299", "6274", "6538", "5744", "5489", "6176", "6797", "6029", "5469", "5188", "5398", "6601", "6227", "5949", "5578", "5980", "5990", "7194", "6398", "6050", "5841", "5840", "6532", "6198", "6258", "6303", "5821", "6552", "6265", "5960", "6238", "6320", "6875", "5919", "6717", "5090", "6792", "4796", "7628", "6979", "6746", "6407", "5944", "6854", "5460", "4806", "6366", "6405", "6507", "6048", "5822", "6946", "5866", "6859", "5123", "6026", "5974", "6962", "5304", "6164", "5772", "6056", "6785", "5853", "6114", "6128", "7128", "6034", "4471", "6130", "6660", "6958", "4744", "7135", "6506", "6620", "6221", "6713", "6393", "5354", "5660", "7273", "6526"
 # First manually copy partion files to path_out
# cp *treefile /home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/20_Astral-mix/10

def move(path_exons, gene_select, path_out, done):
    """, move """
    inputs = [path_out]
    outputs = [done]
    options = {'cores': 2, 'memory': "1g", 'walltime': "00:10:00", 'account':"cryptocarya"}

    spec = """

    cd {path_out}
    rm {gene_select}_part.txt.treefile

    cd {path_exons}
    cp {path_exons}{gene_select}_output_taper.fasta.treefile {path_out}
    cd {path_out}
    mv {path_out}{gene_select}_output_taper.fasta.treefile {path_out}{gene_select}_part.txt.treefile
   
    touch {done}

    """.format(gene_select = gene_select, path_out = path_out, path_exons = path_exons, done = done)

    return (inputs, outputs, options, spec)

# ########################################################################################################################
# #####################################---- Astral Tree Search ----#######################################################
# ########################################################################################################################

# def astral2(path_in, path_out, done):
#     """, ASTRAL """
#     inputs = [path_out]
#     outputs = [done]
#     options = {'cores': 2, 'memory': "10g", 'walltime': "00:30:00", 'account':"cryptocarya"}

#     spec = """

#     # Activate phyx
#     source /home/laurakf/miniconda3/etc/profile.d/conda.sh
#     conda activate Phyx

#     cd {path_in}

#     while read name
#     do cat "$name"_part.txt.treefile >> PAFTOL_trees.tre && echo "", >> PAFTOL_trees.tre
#     done

#     /home/laurakf/cryptocarya/Programs/newick-utils-1.6/src/nw_ed PAFTOL_trees.tre 'i & b<=10' o > PAFTOL_trees_BP10.tre

#     java -jar /home/laurakf/cryptocarya/Programs/Astral/astral.5.7.8.jar -i PAFTOL_trees_BP10.tre -t 2 -o PAFTOL_trees_BP10_SpeciesTree_annotQ.tre

#     java -jar /home/laurakf/cryptocarya/Programs/Astral/astral.5.7.8.jar -i PAFTOL_trees_BP10.tre -t 0 -o PAFTOL_trees_BP10_SpeciesTree.tre

#     pxrr -t PAFTOL_trees_BP10_SpeciesTree.tre -g Myri-fragr-PAFTOL, Magn-grand-PAFTOL > PAFTOL_trees_BP10_SpeciesTree_rooted.tre

#     pxrr -t PAFTOL_trees_BP10_SpeciesTree_annotQ.tre -g Myri-fragr-PAFTOL, Magn-grand-PAFTOL > PAFTOL_trees_BP10_SpeciesTree_annotQ_rooted.tre

#     mv *.tre {path_out}

#     touch {done}

#     """.format(path_out = path_out, path_in = path_in, done = done)

#     return (inputs, outputs, options, spec)


# ########################################################################################################################
# ######################################################---- RUN ----#####################################################
# ########################################################################################################################

# #Species 
# sp = ["Alse-petio-PAFTOL", "Athe-mosch-PAFTOL", "Beil-pendu-PAFTOL", "Beil-tsang-PAFTOL", "Cary-tonki-PAFTOL", "Caly-flori-PAFTOL", "Cass-filif-PAFTOL", "Cinn-camph-PAFTOL", "Cryp-alba-PAFTOL", "Deha-haina-PAFTOL", "Endi-macro-PAFTOL", "Gomo-keule-PAFTOL", "Hern-nymph-PAFTOL", "Idio-austr-PAFTOL", "Laur-nobil-PAFTOL", "Mach-salic-PAFTOL", "Magn-grand-PAFTOL", "Mezi-ita-uba-PAFTOL", "Moll-gilgi-PAFTOL", "Moni-rotun-PAFTOL", "Myri-fragr-PAFTOL", "Neoc-cauda-PAFTOL", "Noth-umbel-PAFTOL", "Pers-borbo-PAFTOL", "Peum-boldu-PAFTOL", "Phoe-lance-PAFTOL", "Sipa-guian-PAFTOL", "Spar-botoc-PAFTOL", "Synd-chine-PAFTOL", "Tamb-ficus-PAFTOL"] 


# for i in range(len(sp)):
#     #### Running fastqc on raw data
#     gwf.target_from_template('fastqc_raw_'+str(i), fastqc_raw(name = sp[i],
#                                                         path_in= "/home/laurakf/cryptocarya/RawData/PAFTOL/",
#                                                         path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/01_FastQC/",
#                                                         done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/01_FastQC/done/"+sp[i]))

# #### Running multiqc on raw data
# gwf.target_from_template('multiqc_raw', multiqc_raw(path_in= "/home/laurakf/cryptocarya/Workflow/PAFTOL/01_FastQC/",
#                                                     path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/02_MultiQC/",
#                                                     done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/02_MultiQC/done/multiqc_raw"))

# sp = ["Alse-petio-PAFTOL", "Athe-mosch-PAFTOL", "Beil-pendu-PAFTOL", "Beil-tsang-PAFTOL", "Cary-tonki-PAFTOL", "Caly-flori-PAFTOL", "Cass-filif-PAFTOL", "Cinn-camph-PAFTOL", "Cryp-alba-PAFTOL", "Deha-haina-PAFTOL", "Endi-macro-PAFTOL", "Gomo-keule-PAFTOL", "Hern-nymph-PAFTOL", "Idio-austr-PAFTOL", "Laur-nobil-PAFTOL", "Mach-salic-PAFTOL", "Magn-grand-PAFTOL", "Mezi-ita-uba-PAFTOL", "Moll-gilgi-PAFTOL", "Moni-rotun-PAFTOL", "Myri-fragr-PAFTOL", "Neoc-cauda-PAFTOL", "Noth-umbel-PAFTOL", "Pers-borbo-PAFTOL", "Peum-boldu-PAFTOL", "Phoe-lance-PAFTOL", "Sipa-guian-PAFTOL", "Spar-botoc-PAFTOL", "Synd-chine-PAFTOL", "Tamb-ficus-PAFTOL"] 

# for i in range(len(sp)):
#     #### Running Trimmomatic (slidingwindow)
#     gwf.target_from_template('trimmomatic_'+str(i), trimmomatic(name = sp[i],
#                                                         path_in= "/home/laurakf/cryptocarya/RawData/PAFTOL/", 
#                                                         path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/03_Trimmomatic/slidingwindow/",
#                                                         done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/03_Trimmomatic/slidingwindow/done/"+sp[i]))

#     #     #### Running Trimmomatic (maxinfo)
#     # gwf.target_from_template('trimmomatic_maxinfo_'+sp[i], trimmomatic_maxinfo(name = sp[i],
#     #                                                     path_in= "/home/laurakf/cryptocarya/Workflow/PAFTOL/01_FastQC/data/", 
#     #                                                     path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/03_Trimmomatic/maxinfo/",
#     #                                                     done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/03_Trimmomatic/maxinfo/done/"+sp[i]))

    # #### Running fastqc on the trimmed data (slidingwindow)
    # gwf.target_from_template('fastqc_trimmed_'+str(i), fastqc_trimmed(name = sp[i],
    #                                                     path_in= "/home/laurakf/cryptocarya/Workflow/PAFTOL/03_Trimmomatic/slidingwindow/secapr_postrim/",
    #                                                     path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/04_FastQC/slidingwindow/",
    #                                                     done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/04_FastQC/slidingwindow/done/"+sp[i]))  

#     # #### Running fastqc on the trimmed data (maxinfo)
#     # gwf.target_from_template('fastqc_trimmed_maxinfo_'+sp[i], fastqc_trimmed_maxinfo(name = sp[i],
#     #                                                     path_in= "/home/laurakf/cryptocarya/Workflow/PAFTOL/03_Trimmomatic/maxinfo/secapr_postrim/",
#     #                                                     path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/04_FastQC/maxinfo/",
#     #                                                     done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/04_FastQC/maxinfo/done/"+sp[i]))                                                  


# #### Running multiqc on trimmed data (slidingwindow)
# gwf.target_from_template('multiqc_trimmed_slidingwindow', multiqc_trimmed(path_in= "/home/laurakf/cryptocarya/Workflow/PAFTOL/04_FastQC/slidingwindow/",
#                                                     path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/05_MultiQC/slidingwindow/",
#                                                     done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/05_MultiQC/slidingwindow/done/multiqc_trimmed"))

# #### Running multiqc on trimmed data (maxinfo)
# gwf.target_from_template('multiqc_trimmed_maxinfo', multiqc_trimmed_maxinfo(path_in= "/home/laurakf/cryptocarya/Workflow/PAFTOL/04_FastQC/maxinfo/",
#                                                     path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/05_MultiQC/maxinfo/",
#                                                     done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/05_MultiQC/maxinfo/done/multiqc_trimmed"))


# sp = ["Alse-petio-PAFTOL", "Athe-mosch-PAFTOL", "Beil-pendu-PAFTOL", "Beil-tsang-PAFTOL", "Cary-tonki-PAFTOL", "Caly-flori-PAFTOL", "Cass-filif-PAFTOL", "Cinn-camph-PAFTOL", "Cryp-alba-PAFTOL", "Deha-haina-PAFTOL", "Endi-macro-PAFTOL", "Gomo-keule-PAFTOL", "Hern-nymph-PAFTOL", "Idio-austr-PAFTOL", "Laur-nobil-PAFTOL", "Mach-salic-PAFTOL", "Magn-grand-PAFTOL", "Mezi-ita-uba-PAFTOL", "Moll-gilgi-PAFTOL", "Moni-rotun-PAFTOL", "Myri-fragr-PAFTOL", "Neoc-cauda-PAFTOL", "Noth-umbel-PAFTOL", "Pers-borbo-PAFTOL", "Phoe-lance-PAFTOL", "Sipa-guian-PAFTOL", "Spar-botoc-PAFTOL", "Tamb-ficus-PAFTOL"] 
# # Taken Synd-chine-PAFTOL out = too large. Taken Peum-boldu-PAFTOL out. They do not seem to work. 

# for i in range(len(sp)):
    
#     #### Running Hybpiper
#     gwf.target_from_template('Hybpiper_'+str(i), hybpiper(name = sp[i],
#                                                         p1 = "_1P.fastq",
#                                                         p2 = "_2P.fastq",
#                                                         un = "_UN.fastq",
#                                                         path_out= "/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/",
#                                                         path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/03_Trimmomatic/slidingwindow/",
#                                                         done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/done/HybPiper/"+sp[i]))

# #### Getting stats and heatmap
# gwf.target_from_template('stats', stats(path_out= "/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/Stats_Heatmap/",
#                                                 path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/",
#                                                 in_done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/done/HybPiper/"+sp[i],
#                                                 name = sp[i],
#                                                 done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/done/Stats/"+sp[i]))

                                               
# #### Paralogs
# gwf.target_from_template('Paralogs', paralogs(name = sp[i],
#                                                       path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/",
#                                                       path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/Paralogs/", ,
#                                                       done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/done/Paralogs/done",
#                                                       in_done="/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/done/HybPiper/"+sp[i]))

# sp = ["Alse-petio-PAFTOL", "Athe-mosch-PAFTOL", "Beil-pendu-PAFTOL", "Beil-tsang-PAFTOL", "Cary-tonki-PAFTOL", "Caly-flori-PAFTOL", "Cass-filif-PAFTOL", "Cinn-camph-PAFTOL", "Cryp-alba-PAFTOL", "Deha-haina-PAFTOL", "Endi-macro-PAFTOL", "Gomo-keule-PAFTOL", "Hern-nymph-PAFTOL", "Idio-austr-PAFTOL", "Laur-nobil-PAFTOL", "Mach-salic-PAFTOL", "Magn-grand-PAFTOL", "Mezi-ita-uba-PAFTOL", "Moll-gilgi-PAFTOL", "Moni-rotun-PAFTOL", "Myri-fragr-PAFTOL", "Neoc-cauda-PAFTOL", "Noth-umbel-PAFTOL", "Pers-borbo-PAFTOL", "Phoe-lance-PAFTOL", "Sipa-guian-PAFTOL", "Spar-botoc-PAFTOL", "Tamb-ficus-PAFTOL"] 
# # Taken Synd-chine-PAFTOL out = too large. Taken Peum-boldu-PAFTOL out. They do not seem to work. 

# for i in range(len(sp)):

#     #### Coverage
#     gwf.target_from_template('Coverage_'+str(i), coverage(name = sp[i],
#                                                         path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/",
#                                                         all_bam = "_all.bam",
#                                                         all_sorted_bam ="_all_sorted.bam",
#                                                         all_sorted_bam_bai="_all_sorted.bam.bai",
#                                                         bam =".bam",
#                                                         cov=".cov",
#                                                         fasta = ".fasta",
#                                                         fasta_amb = ".fasta.amb",
#                                                         fasta_ann = ".fasta.ann",
#                                                         fasta_bwt = ".fasta.bwt",
#                                                         fasta_pac = ".fasta.pac",
#                                                         fasta_sa = ".fasta.sa",
#                                                         trimmed_fasta = "_trimmed.fasta",
#                                                         up_bam = "_up.bam",
#                                                         path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/07_Coverage/",
#                                                         done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/07_Coverage/done/Coverage/"+sp[i],
#                                                         dir_wrk = "/home/laurakf/cryptocarya/Workflow/PAFTOL/06_HybPiper/",
#                                                         dir_in ="/home/laurakf/cryptocarya/Workflow/PAFTOL/03_Trimmomatic/slidingwindow/", #Folder with clean reads + unpaired
#                                                         dir_out ="/home/laurakf/cryptocarya/Workflow/PAFTOL/07_Coverage/")) # folder with coverage

# #### Retrieve sequences and sort into files with gene names
# gwf.target_from_template('retrieve', retrieve(path_in ="/home/laurakf/cryptocarya/Workflow/PAFTOL/07_Coverage/", 
#                                               done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/08_Retrieve/done/retrieve"))


# ### Kan jeg vente med at fjerne de gener der måske er paraloge indtil dette step, hvor jeg vælger hvilke gener der skal inkluderes? ### JA
# ### Skal jeg fjerne ITS og kloroplast regionerne i den videre analyse? JA

# gene = ["4471", "4527", "4691", "4724", "4744", "4757", "4793", "4796", "4802", "4806", "4848", "4889", "4890", "4893", "4932", "4942", "4951", "4954", "4992", "5018", "5032", "5034", "5038", "5090", "5104", "5116", "5123", "5131", "5138", "5162", "5163", "5168", "5177", "5188", "5200", "5206", "5220", "5257", "5264", "5271", "5273", "5280", "5296", "5299", "5304", "5318", "5326", "5328", "5333", "5335", "5339", "5343", "5348", "5354", "5366", "5398", "5404", "5406", "5421", "5426", "5427", "5428", "5449", "5454", "5460", "5463", "5464", "5469", "5477", "5489", "5502", "5513", "5528", "5531", "5536", "5551", "5554", "5562", "5578", "5594", "5596", "5599", "5614", "5620", "5634", "5639", "5644", "5656", "5660", "5664", "5670", "5699", "5702", "5703", "5716", "5721", "5744", "5770", "5772", "5791", "5802", "5815", "5816", "5821", "5822", "5840", "5841", "5842", "5843", "5849", "5853", "5857", "5858", "5865", "5866", "5870", "5893", "5894", "5899", "5910", "5913", "5918", "5919", "5921", "5926", "5933", "5942", "5943", "5944", "5945", "5949", "5950", "5960", "5968", "5974", "5977", "5980", "5981", "5990", "6000", "6003", "6004", "6016", "6026", "6029", "6034", "6036", "6038", "6041", "6048", "6050", "6051", "6056", "6064", "6068", "6072", "6098", "6114", "6119", "6128", "6130", "6139", "6148", "6164", "6175", "6176", "6198", "6216", "6221", "6226", "6227", "6238", "6258", "6265", "6274", "6282", "6284", "6295", "6298", "6299", "6303", "6318", "6320", "6363", "6366", "6373", "6376", "6378", "6383", "6384", "6389", "6393", "6398", "6404", "6405", "6406", "6407", "6412", "6420", "6432", "6439", "6447", "6450", "6454", "6457", "6458", "6459", "6460", "6462", "6483", "6487", "6488", "6492", "6494", "6496", "6500", "6506", "6507", "6526", "6527", "6528", "6532", "6533", "6538", "6540", "6544", "6550", "6552", "6559", "6563", "6570", "6572", "6601", "6620", "6631", "6636", "6639", "6641", "6649", "6652", "6660", "6667", "6685", "6689", "6713", "6717", "6732", "6733", "6738", "6746", "6779", "6782", "6785", "6792", "6797", "6825", "6848", "6854", "6859", "6860", "6865", "6875", "6882", "6883", "6909", "6913", "6914", "6924", "6933", "6946", "6947", "6954", "6958", "6961", "6962", "6968", "6978", "6979", "6992", "7021", "7024", "7029", "7067", "7111", "7128", "7135", "7136", "7141", "7174", "7194", "7241", "7273", "7279", "7313", "7324", "7325", "7331", "7333", "7336", "7363", "7367", "7371", "7572", "7577", "7583", "7602", "7628"]
# genes = ["4471.FNA", "4527.FNA", "4691.FNA", "4724.FNA", "4744.FNA", "4757.FNA", "4793.FNA", "4796.FNA", "4802.FNA", "4806.FNA", "4848.FNA", "4889.FNA", "4890.FNA", "4893.FNA", "4932.FNA", "4942.FNA", "4951.FNA", "4954.FNA", "4992.FNA", "5018.FNA", "5032.FNA", "5034.FNA", "5038.FNA", "5090.FNA", "5104.FNA", "5116.FNA", "5123.FNA", "5131.FNA", "5138.FNA", "5162.FNA", "5163.FNA", "5168.FNA", "5177.FNA", "5188.FNA", "5200.FNA", "5206.FNA", "5220.FNA", "5257.FNA", "5264.FNA", "5271.FNA", "5273.FNA", "5280.FNA", "5296.FNA", "5299.FNA", "5304.FNA", "5318.FNA", "5326.FNA", "5328.FNA", "5333.FNA", "5335.FNA", "5339.FNA", "5343.FNA", "5348.FNA", "5354.FNA", "5366.FNA", "5398.FNA", "5404.FNA", "5406.FNA", "5421.FNA", "5426.FNA", "5427.FNA", "5428.FNA", "5449.FNA", "5454.FNA", "5460.FNA", "5463.FNA", "5464.FNA", "5469.FNA", "5477.FNA", "5489.FNA", "5502.FNA", "5513.FNA", "5528.FNA", "5531.FNA", "5536.FNA", "5551.FNA", "5554.FNA", "5562.FNA", "5578.FNA", "5594.FNA", "5596.FNA", "5599.FNA", "5614.FNA", "5620.FNA", "5634.FNA", "5639.FNA", "5644.FNA", "5656.FNA", "5660.FNA", "5664.FNA", "5670.FNA", "5699.FNA", "5702.FNA", "5703.FNA", "5716.FNA", "5721.FNA", "5744.FNA", "5770.FNA", "5772.FNA", "5791.FNA", "5802.FNA", "5815.FNA", "5816.FNA", "5821.FNA", "5822.FNA", "5840.FNA", "5841.FNA", "5842.FNA", "5843.FNA", "5849.FNA", "5853.FNA", "5857.FNA", "5858.FNA", "5865.FNA", "5866.FNA", "5870.FNA", "5893.FNA", "5894.FNA", "5899.FNA", "5910.FNA", "5913.FNA", "5918.FNA", "5919.FNA", "5921.FNA", "5926.FNA", "5933.FNA", "5942.FNA", "5943.FNA", "5944.FNA", "5945.FNA", "5949.FNA", "5950.FNA", "5960.FNA", "5968.FNA", "5974.FNA", "5977.FNA", "5980.FNA", "5981.FNA", "5990.FNA", "6000.FNA", "6003.FNA", "6004.FNA", "6016.FNA", "6026.FNA", "6029.FNA", "6034.FNA", "6036.FNA", "6038.FNA", "6041.FNA", "6048.FNA", "6050.FNA", "6051.FNA", "6056.FNA", "6064.FNA", "6068.FNA", "6072.FNA", "6098.FNA", "6114.FNA", "6119.FNA", "6128.FNA", "6130.FNA", "6139.FNA", "6148.FNA", "6164.FNA", "6175.FNA", "6176.FNA", "6198.FNA", "6216.FNA", "6221.FNA", "6226.FNA", "6227.FNA", "6238.FNA", "6258.FNA", "6265.FNA", "6274.FNA", "6282.FNA", "6284.FNA", "6295.FNA", "6298.FNA", "6299.FNA", "6303.FNA", "6318.FNA", "6320.FNA", "6363.FNA", "6366.FNA", "6373.FNA", "6376.FNA", "6378.FNA", "6383.FNA", "6384.FNA", "6389.FNA", "6393.FNA", "6398.FNA", "6404.FNA", "6405.FNA", "6406.FNA", "6407.FNA", "6412.FNA", "6420.FNA", "6432.FNA", "6439.FNA", "6447.FNA", "6450.FNA", "6454.FNA", "6457.FNA", "6458.FNA", "6459.FNA", "6460.FNA", "6462.FNA", "6483.FNA", "6487.FNA", "6488.FNA", "6492.FNA", "6494.FNA", "6496.FNA", "6500.FNA", "6506.FNA", "6507.FNA", "6526.FNA", "6527.FNA", "6528.FNA", "6532.FNA", "6533.FNA", "6538.FNA", "6540.FNA", "6544.FNA", "6550.FNA", "6552.FNA", "6559.FNA", "6563.FNA", "6570.FNA", "6572.FNA", "6601.FNA", "6620.FNA", "6631.FNA", "6636.FNA", "6639.FNA", "6641.FNA", "6649.FNA", "6652.FNA", "6660.FNA", "6667.FNA", "6685.FNA", "6689.FNA", "6713.FNA", "6717.FNA", "6732.FNA", "6733.FNA", "6738.FNA", "6746.FNA", "6779.FNA", "6782.FNA", "6785.FNA", "6792.FNA", "6797.FNA", "6825.FNA", "6848.FNA", "6854.FNA", "6859.FNA", "6860.FNA", "6865.FNA", "6875.FNA", "6882.FNA", "6883.FNA", "6909.FNA", "6913.FNA", "6914.FNA", "6924.FNA", "6933.FNA", "6946.FNA", "6947.FNA", "6954.FNA", "6958.FNA", "6961.FNA", "6962.FNA", "6968.FNA", "6978.FNA", "6979.FNA", "6992.FNA", "7021.FNA", "7024.FNA", "7029.FNA", "7067.FNA", "7111.FNA", "7128.FNA", "7135.FNA", "7136.FNA", "7141.FNA", "7174.FNA", "7194.FNA", "7241.FNA", "7273.FNA", "7279.FNA", "7313.FNA", "7324.FNA", "7325.FNA", "7331.FNA", "7333.FNA", "7336.FNA", "7363.FNA", "7367.FNA", "7371.FNA", "7572.FNA", "7577.FNA", "7583.FNA", "7602.FNA", "7628.FNA"]


# #### MAFFT
# for i in range(len(genes)):
#     gwf.target_from_template('Mafft_'+genes[i], mafft(genes = genes[i],
#                                                         path_out= "/home/laurakf/cryptocarya/Workflow/PAFTOL/09_Mafft/",
#                                                         path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/08_Retrieve/",
#                                                         gene = gene[i],
#                                                         done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/09_Mafft/done/"+genes[i]))


# #### Trimal according to a pre-defined gt values
# for i in range(len(gene)):
#    gwf.target_from_template('gt_trimming_'+gene[i], gt_trimming(gene = gene[i],
#                                                         path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/09_Mafft/",
#                                                         path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/10_Trimal/",
#                                                         done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/10_Trimal/done/"+gene[i]))

# ### Generating AMAS statistics for raw_alignments
# gwf.target_from_template('amas_raw', amas_raw(path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/10_Trimal",
#                                         in_done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/10_Trimal/done/",
#                                         done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/10_Trimal/done/AMAS_raw/raw"))

# cut_off = ["0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9"]

# #### Generating AMAS statistics for gt_alignments
# for i in range(len(cut_off)):
#     gwf.target_from_template('amas_gt_'+cut_off[i], amas_gt(path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/10_Trimal/",
#                                                 cut_off = cut_off[i],
#                                                 in_done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/10_Trimal/done/"+gene[i],
#                                                 done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/10_Trimal/done/AMAS_gt/"+cut_off[i]))


# #### Optrimal
# gwf.target_from_template('optrim', optrim(path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/10_Trimal/",
#                                                          done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/11_Optrimal/done/optrimal/optrimal", 
#                                                          path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/11_Optrimal/"))


# #### Move files from trimal to optrimal folder
# gwf.target_from_template('move', move(path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/10_Trimal/",
#                                                 in_done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/11_Optrimal/done/optrimal/optrimal",
#                                                 done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/11_Optrimal/done/move/move", 
#                                                 path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/11_Optrimal/"))

# Removing genes that were not included in "optimal_final_results", - 
# gene = ["4471", "4527", "4691", "4724", "4744", "4757", "4793", "4796", "4802", "4806", "4848", "4889", "4890", "4893", "4932", "4942", "4951", "4954", "4992", "5018", "5032", "5034", "5038", "5090", "5104", "5116", "5123", "5131", "5138", "5162", "5163", "5168", "5177", "5188", "5200", "5206", "5220", "5257", "5264", "5271", "5273", "5280", "5296", "5299", "5304", "5318", "5326", "5328", "5333", "5335", "5339", "5343", "5348", "5354", "5366", "5398", "5404", "5406", "5421", "5426", "5427", "5428", "5449", "5454", "5460", "5463", "5464", "5469", "5477", "5489", "5502", "5513", "5528", "5531", "5536", "5551", "5554", "5562", "5578", "5594", "5596", "5614", "5620", "5634", "5639", "5644", "5656", "5660", "5664", "5670", "5699", "5702", "5703", "5716", "5721", "5744", "5770", "5772", "5791", "5802", "5815", "5816", "5821", "5822", "5840", "5841", "5842", "5843", "5849", "5853", "5857", "5858", "5865", "5866", "5870", "5893", "5894", "5899", "5910", "5913", "5918", "5919", "5926", "5933", "5942", "5944", "5945", "5949", "5950", "5960", "5968", "5974", "5977", "5980", "5981", "5990", "6000", "6003", "6004", "6016", "6026", "6029", "6034", "6036", "6038", "6048", "6050", "6051", "6056", "6064", "6068", "6072", "6098", "6114", "6119", "6128", "6130", "6139", "6148", "6164", "6175", "6176", "6198", "6216", "6221", "6226", "6227", "6238", "6258", "6265", "6274", "6282", "6284", "6295", "6298", "6299", "6303", "6318", "6320", "6363", "6366", "6373", "6376", "6378", "6383", "6384", "6389", "6393", "6398", "6404", "6405", "6406", "6407", "6412", "6420", "6432", "6439", "6447", "6450", "6454", "6457", "6458", "6459", "6460", "6462", "6483", "6487", "6488", "6492", "6494", "6496", "6500", "6506", "6507", "6526", "6527", "6528", "6532", "6533", "6538", "6540", "6544", "6550", "6552", "6559", "6563", "6570", "6572", "6601", "6620", "6631", "6636", "6639", "6641", "6649", "6652", "6660", "6667", "6685", "6689", "6713", "6717", "6732", "6733", "6738", "6746", "6779", "6782", "6785", "6792", "6797", "6825", "6848", "6854", "6859", "6860", "6865", "6875", "6882", "6883", "6909", "6913", "6914", "6924", "6933", "6946", "6947", "6954", "6958", "6961", "6962", "6968", "6978", "6979", "6992", "7021", "7024", "7029", "7067", "7111", "7128", "7135", "7136", "7141", "7174", "7194", "7241", "7273", "7279", "7313", "7324", "7325", "7331", "7333", "7336", "7363", "7367", "7371", "7572", "7577", "7583", "7602", "7628"]                                            
# Running CIAlign on the trimmed_fasta - excluded some genes earlier with paralog warnings.
# for i in range(0, len(gene)):
#     # gwf.target_from_template('Cialign'+gene[i], cialign1(gene = gene[i],
    #                                           path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/11_Optrimal/optimal_final_results/",
    #                                           path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/12_CIAlign/",
    #                                           done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/12_CIAlign/done/"+gene[i]))

# ## Running TAPER after CIALIGN
#     gwf.target_from_template('Taper_'+gene[i], taper(gene = gene[i],
#                                                     path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/12_CIAlign/",
#                                                     path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL/13_Taper/",
#                                                     done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/13_Taper/done/"+gene[i]))

    # #### Running Exon_mapper
    # gwf.target_from_template('Exon_map_'+gene[i], exon_map(gene = gene[i],
    #                                                     path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/13_Taper/",
    #                                                     path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/17_mapping/",
    #                                                     done = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/17_mapping/done/"+gene[i]))

   # ### Creating the partition files for each gene
   #  gwf.target_from_template('Partition_'+gene[i], partitioner(gene = gene[i],
   #                                                      path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/17_mapping/",
   #                                                      path_out= "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/18_partitions/",
   #                                                      done = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/18_partitions/done/"+gene[i]))


# gene = ["4471", "4527", "4691", "4724", "4744", "4757", "4793", "4796", "4802", "4806", "4848", "4889", "4890", "4893", "4932", "4942", "4951", "4954", "4992", "5018", "5032", "5034", "5038", "5090", "5104", "5116", "5123", "5131", "5138", "5162", "5163", "5168", "5177", "5188", "5200", "5206", "5220", "5257", "5264", "5271", "5273", "5280", "5296", "5299", "5304", "5318", "5326", "5328", "5333", "5335", "5339", "5343", "5348", "5354", "5366", "5398", "5404", "5406", "5421", "5426", "5427", "5428", "5449", "5454", "5460", "5463", "5464", "5469", "5477", "5489", "5502", "5513", "5528", "5531", "5536", "5551", "5554", "5562", "5578", "5594", "5596", "5614", "5620", "5634", "5639", "5644", "5656", "5660", "5664", "5670", "5699", "5702", "5703", "5716", "5721", "5744", "5770", "5772", "5791", "5802", "5815", "5816", "5821", "5822", "5840", "5841", "5842", "5843", "5849", "5853", "5857", "5858", "5865", "5866", "5870", "5893", "5894", "5899", "5910", "5913", "5918", "5919", "5926", "5933", "5942", "5944", "5945", "5949", "5950", "5960", "5968", "5974", "5977", "5980", "5981", "5990", "6000", "6003", "6004", "6016", "6026", "6029", "6034", "6036", "6038", "6048", "6050", "6051", "6056", "6064", "6068", "6072", "6098", "6114", "6119", "6128", "6130", "6139", "6148", "6164", "6175", "6176", "6198", "6216", "6221", "6226", "6227", "6238", "6258", "6265", "6274", "6282", "6284", "6295", "6298", "6299", "6303", "6318", "6320", "6363", "6366", "6373", "6376", "6378", "6383", "6384", "6389", "6393", "6398", "6404", "6405", "6406", "6407", "6412", "6420", "6432", "6439", "6447", "6450", "6454", "6457", "6458", "6459", "6460", "6462", "6483", "6487", "6488", "6492", "6494", "6496", "6500", "6506", "6507", "6526", "6527", "6528", "6532", "6533", "6538", "6540", "6544", "6550", "6552", "6559", "6563", "6570", "6572", "6601", "6620", "6631", "6636", "6639", "6641", "6649", "6652", "6660", "6667", "6685", "6689", "6713", "6717", "6732", "6733", "6738", "6746", "6779", "6782", "6785", "6792", "6797", "6825", "6848", "6854", "6859", "6860", "6865", "6875", "6882", "6883", "6909", "6913", "6914", "6924", "6933", "6946", "6947", "6954", "6958", "6961", "6962", "6968", "6978", "6979", "6992", "7021", "7024", "7029", "7067", "7111", "7128", "7135", "7136", "7141", "7174", "7194", "7241", "7273", "7279", "7313", "7324", "7325", "7331", "7333", "7336", "7363", "7367", "7371", "7572", "7577", "7583", "7602", "7628"]                                            


# #Running IQTREE for files trimmed with trimal and CIAlign                                             
# for i in range(0, len(gene)):
#    gwf.target_from_template('Iqtree_'+gene[i], iq_tree(gene = gene[i],
#                                                     path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/14_IQtree/",
#                                                     done = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/14_IQtree/done/"+gene[i],
#                                                     path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/18_partitions/"))  


# # Running ASTRAL 
# gwf.target_from_template('astral_tapper', astral(path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL/15_Raxml/",
#                                                     done = "/home/laurakf/cryptocarya/Workflow/PAFTOL/16_Astral/done/Astral"))

# gene = ["4471", "4527", "4691", "4724", "4744", "4757", "4793", "4796", "4802", "4806", "4848", "4889", "4890", "4932", "4942", "4951", "4954", "4992", "5018", "5032", "5034", "5038", "5090", "5104", "5116", "5123", "5131", "5138", "5162", "5163", "5168", "5177", "5188", "5200", "5206", "5220", "5257", "5264", "5271", "5273", "5280", "5296", "5299", "5304", "5318", "5326", "5328", "5333", "5335", "5339", "5343", "5348", "5354", "5366", "5398", "5404", "5406", "5421", "5426", "5427", "5428", "5449", "5454", "5460", "5463", "5464", "5469", "5477", "5489", "5502", "5513", "5528", "5531", "5536", "5551", "5554", "5562", "5578", "5594", "5596", "5614", "5620", "5634", "5639", "5644", "5656", "5660", "5664", "5670", "5699", "5702", "5703", "5716", "5721", "5744", "5770", "5772", "5791", "5802", "5815", "5816", "5821", "5822", "5840", "5841", "5842", "5843", "5849", "5853", "5857", "5858", "5865", "5866", "5870", "5893", "5894", "5899", "5910", "5913", "5918", "5919", "5926", "5933", "5942", "5944", "5945", "5949", "5950", "5960", "5968", "5974", "5977", "5980", "5981", "5990", "6000", "6003", "6004", "6016", "6026", "6029", "6034", "6036", "6038", "6048", "6050", "6051", "6056", "6064", "6068", "6072", "6098", "6114", "6119", "6128", "6130", "6139", "6148", "6164", "6175", "6176", "6198", "6216", "6221", "6226", "6227", "6238", "6258", "6265", "6274", "6282", "6284", "6295", "6298", "6299", "6303", "6318", "6320", "6363", "6366", "6373", "6376", "6378", "6383", "6384", "6389", "6393", "6398", "6404", "6405", "6407", "6412", "6420", "6432", "6439", "6447", "6450", "6454", "6457", "6458", "6459", "6460", "6462", "6483", "6488", "6492", "6494", "6496", "6500", "6506", "6507", "6526", "6527", "6528", "6532", "6533", "6538", "6540", "6544", "6550", "6552", "6559", "6563", "6570", "6572", "6601", "6620", "6631", "6636", "6639", "6641", "6649", "6652", "6660", "6667", "6685", "6689", "6713", "6717", "6732", "6733", "6738", "6746", "6779", "6782", "6785", "6792", "6797", "6825", "6848", "6854", "6859", "6860", "6865", "6875", "6882", "6883", "6909", "6913", "6914", "6924", "6933", "6946", "6947", "6954", "6958", "6961", "6962", "6968", "6978", "6979", "6992", "7021", "7024", "7029", "7067", "7111", "7128", "7135", "7136", "7141", "7174", "7194", "7241", "7273", "7279", "7313", "7324", "7325", "7331", "7333", "7336", "7363", "7367", "7371", "7572", "7577", "7583", "7602", "7628"]                                            
# # Removed 3 genes not present in rooted folder (4893, 6487, 7628)

# for i in range(0, len(gene)):
#     #Running Root gene trees
#     gwf.target_from_template('Root_'+gene[i], root_genetrees(gene = gene[i],
#                                                         path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/19_sortadate/rooted_geneTrees/",
#                                                         path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/19_sortadate/gene_trees/"))

# # Running SortaDate on the Astral tree using the gene trees
# gwf.target_from_template('Sorta_date_partition', sorta_date(path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/19_sortadate/rooted_geneTrees/",
#                                                         path_out ="/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/19_sortadate/output/",
#                                                         astral_tree="/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/16_Astral/PAFTOL_trees_BP10_SpeciesTree_rooted2.tre",
#                                                         done = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/19_sortadate/done/sorted"))

# # Exon genes from Taper output
gene = ["4471", "4527", "4691", "4724", "4744", "4757", "4793", "4796", "4802", "4806", "4848", "4889", "4890", "4893", "4932", "4942", "4951", "4954", "4992", "5018", "5032", "5034", "5038", "5090", "5104", "5116", "5123", "5131", "5138", "5162", "5163", "5168", "5177", "5188", "5200", "5206", "5220", "5257", "5264", "5271", "5273", "5280", "5296", "5299", "5304", "5318", "5326", "5328", "5333", "5335", "5339", "5343", "5348", "5354", "5366", "5398", "5404", "5406", "5421", "5426", "5427", "5428", "5449", "5454", "5460", "5463", "5464", "5469", "5477", "5489", "5502", "5513", "5528", "5531", "5536", "5551", "5554", "5562", "5578", "5594", "5596", "5614", "5620", "5634", "5639", "5644", "5656", "5660", "5664", "5670", "5699", "5702", "5703", "5716", "5721", "5744", "5770", "5772", "5791", "5802", "5815", "5816", "5821", "5822", "5840", "5841", "5842", "5843", "5849", "5853", "5857", "5858", "5865", "5866", "5870", "5893", "5894", "5899", "5910", "5913", "5918", "5919", "5926", "5933", "5942", "5944", "5945", "5949", "5950", "5960", "5968", "5974", "5977", "5980", "5981", "5990", "6000", "6003", "6004", "6016", "6026", "6029", "6034", "6036", "6038", "6048", "6050", "6051", "6056", "6064", "6068", "6072", "6098", "6114", "6119", "6128", "6130", "6139", "6148", "6164", "6175", "6176", "6198", "6216", "6221", "6226", "6227", "6238", "6258", "6265", "6274", "6282", "6284", "6295", "6298", "6299", "6303", "6318", "6320", "6363", "6366", "6373", "6376", "6378", "6383", "6384", "6389", "6393", "6398", "6404", "6405", "6406", "6407", "6412", "6420", "6432", "6439", "6447", "6450", "6454", "6457", "6458", "6459", "6460", "6462", "6483", "6487", "6488", "6492", "6494", "6496", "6500", "6506", "6507", "6526", "6527", "6528", "6532", "6533", "6538", "6540", "6544", "6550", "6552", "6559", "6563", "6570", "6572", "6601", "6620", "6631", "6636", "6639", "6641", "6649", "6652", "6660", "6667", "6685", "6689", "6713", "6717", "6732", "6733", "6738", "6746", "6779", "6782", "6785", "6792", "6797", "6825", "6848", "6854", "6859", "6860", "6865", "6875", "6882", "6883", "6909", "6913", "6914", "6924", "6933", "6946", "6947", "6954", "6958", "6961", "6962", "6968", "6978", "6979", "6992", "7021", "7024", "7029", "7067", "7111", "7128", "7135", "7136", "7141", "7174", "7194", "7241", "7273", "7279", "7313", "7324", "7325", "7331", "7333", "7336", "7363", "7367", "7371", "7572", "7577", "7583", "7602", "7628"]                                            

# #Running IQTREE for exon files trimmed with trimal and CIAlign (to get same format as supercontig partitioner)                                            
# for i in range(0, len(gene)):
#    gwf.target_from_template('Iqtree_'+gene[i], iq_tree_exon(gene = gene[i],
#                                                     path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL-exons/17_IQtree_geneTrees",
#                                                     done = "/home/laurakf/cryptocarya/Workflow/PAFTOL-exons/17_IQtree_geneTrees/done/"+gene[i],
#                                                     path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL-exons/13_Taper/"))  

gene_select = ["4848", "6175", "6282", "6038", "6412", "6483", "6667", "5464", "5264", "5664", "6860", "6882", "5865", "5104", "5670", "5220", "6733", "6636", "5703", "5551", "5950", "6389", "6559", "5200", "7241", "5644", "6572", "5018", "5721", "7331", "6454", "5366", "6685", "6439", "6825", "5257", "5528", "7367", "5296", "5206", "4890", "6000", "6298", "6492", "5918", "5843", "5273", "5802", "7325", "7333", "6299", "6378", "5536", "5513", "4889", "6527", "5933", "6004", "5968", "6550", "6284", "6383", "5531", "6148", "6226", "5299", "6274", "6538", "5744", "5489", "6176", "6797", "6029", "5469", "5188", "5398", "6601", "6227", "5949", "5578", "5980", "5990", "7194", "6398", "6050", "5841", "5840", "6532", "6198", "6258", "6303", "5821", "6552", "6265", "5960", "6238", "6320", "6875", "5919", "6717", "5090", "6792", "4796", "7628", "6979", "6746", "6407", "5944", "6854", "5460", "4806", "6366", "6405", "6507", "6048", "5822", "6946", "5866", "6859", "5123", "6026", "5974", "6962", "5304", "6164", "5772", "6056", "6785", "5853", "6114", "6128", "7128", "6034", "4471", "6130", "6660", "6958", "4744", "7135", "6506", "6620", "6221", "6713", "6393", "5354", "5660", "7273", "6526"]
for i in range(0, len(gene_select)):
   gwf.target_from_template('gene_select_'+gene_select[i], move(gene_select = gene_select[i],
                                                    path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/20_Astral-mix/50/",
                                                    done = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/20_Astral-mix/50/done/gene_select/"+gene_select[i],
                                                    path_exons = "/home/laurakf/cryptocarya/Workflow/PAFTOL-exons/17_IQtree_geneTrees/gene_trees/"))  





# gwf.target_from_template('astral_', astral2(path_out = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/20_Astral-mix/10/astral/",
#                                             done = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/20_Astral-mix/10/done/astral/astral",
#                                             path_in = "/home/laurakf/cryptocarya/Workflow/PAFTOL-partitions/20_Astral-mix/10/"))  


