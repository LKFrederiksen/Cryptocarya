#!/bin/bash
#
for f in *.fasta; do (sed -i'.old' -e 's/n/-/g' $f); done


#Getting the gene name as a commandline argument with the flag -g
while getopts g: flag
do
    case "${flag}" in
        g) gene=${OPTARG};;
    esac
done

#Activating trimal_env
    source /home/laurakf/miniconda3/etc/profile.d/conda.sh
    conda activate Trimal

    # create summary tables for all thresholds specified
    while read cutoff_trim
    do
                        # Checking if a directory exists for the cutoff_trim
            if [[ -d "/home/laurakf/cryptocarya/Workflow/Final_tree/10_Trimal/]${cutoff_trim}" ]]
                        then
                        echo "${cutoff_trim} folder exists."
                        else
                                mkdir $cutoff_trim
                        fi

                        #trimming the aligned sequences of a gene with the given cutoff_trim
            for alignment in $gene
            do
              trimal -in ${alignment} -out ${cutoff_trim}/${alignment} -htmlout ${cutoff_trim}/${alignment/.fasta}.html -gt ${cutoff_trim}

                    #check if alignment was trimmed to extinction by trimAl
                    if grep ' 0 bp' ${cutoff_trim}/${alignment}
                    then
                        rm -f ${cutoff_trim}/${alignment}
                    fi
            done

    done < cutoff_trim.txt
