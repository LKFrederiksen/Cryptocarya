#!/usr/bin/python3

import argparse, dendropy, subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--treefile")
parser.add_argument("--gene")

args = parser.parse_args()
treefile = str(args.treefile)
gene = args.gene

tree = dendropy.Tree.get(path=treefile, schema="newick")
tips = str(tree.taxon_namespace)

#Rooting tree on outgroup species. 
# Myri-fragr-PAFTOL, Myristica fragrans
# Magn-grand-PAFTOL, Magnolia grandiflora


if "Myri-fragr-PAFTOL" in tips and "Magn-grand-PAFTOL" in tips:
    cmd = "pxrr -t "+treefile+" -g Myri-fragr-PAFTOL,Magn-grand-PAFTOL -o {gene}_rooted.tre".format(gene=gene)
else:
    if  "Myri-fragr-PAFTOL" in tips: 
        cmd = "pxrr -t "+treefile+" -g Myri-fragr-PAFTOL -o {gene}_rooted.tre".format(gene=gene)
    elif "Magn-grand-PAFTOL" in tips:
        cmd = "pxrr -t "+treefile+" -g Magn-grand-PAFTOL -o {gene}_rooted.tre".format(gene=gene)
    else:
        cmd = "mv "+treefile+" /home/laurakf/cryptocarya/Workflow/Final_tree/18_SortaDate/no_outgroup/"
        raise ValueError("Outgroup not in ", treefile)

subprocess.call(cmd, shell=True)