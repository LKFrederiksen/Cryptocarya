#!/usr/bin/env python3

import os, subprocess, argparse
from Bio import SeqIO

# parameter for ignoring mini-partitions <smoother
parser = argparse.ArgumentParser()
parser.add_argument("--smoother", default=10, help='minimum length of a partition to be accepted (default 10bp)')
parser.add_argument("--gene", help="The gene for which you want to create a partition")
parser.add_argument("--file_ending", help="The file endings of the fasta files you want to create partitions for.")
args = parser.parse_args()
smoother = int(args.smoother)
gene = args.gene
file_ending = args.file_ending


sequences = [] # gather sequences to keep in the final alignment (all but the exons)
# extract aligned exon sequences
for record in SeqIO.parse(gene+file_ending, "fasta"):
	print(record)
	if record.id == "exon1":
		exon1 = record.seq
	elif record.id == "exon2":
		exon2 = record.seq
	else:
		sequences.append(record)
# create binary partition (1 = exon, 0 = intron)
binpart = []
for i in range(len(exon1)):
	# assumes that any alignment pos. where ANY of the two exons has a base is exon
	if  len(exon1) == len(exon2):
		if exon1[i] == "-" and exon2[i] == "-":
			binpart.append(0)
		else:
			binpart.append(1)
	else:
			print("Error", gene ,len(exon1),len(exon2))
			badgenes = open("bad_exons.txt","a")
			badgenes.write(gene)
			badgenes.write(str(len(exon1)))
			badgenes.write(str(len(exon2+"\n")))
			badgenes.close()
			break
			
	# define partitions as ranges
	partitions_intron = []
	partitions_exon = []
	current_part = [1, "NA"]
	current_state = binpart[0]
	for i in range(len(binpart)):
		if binpart[i] != current_state: # when a shift occurs...
			# define search window for smoothing over mini-partitions
			if (i+smoother) < len(binpart): # to avoid index errors
				wndw = binpart[i:(i+smoother)]
			else:
				wndw = binpart[i:]
			if len(set(wndw)) == 1: # only invoke a partition shift if there is no further partition shift within n=smoother basepairs upstream
				current_part[1] = i # set end position of current partition
				if current_state == 0:  # append concluded partition to partitions_intron or partitions_exon, depending on current_state
					partitions_intron.append(current_part)
				else:
					partitions_exon.append(current_part)
				current_part = [i+1, "NA"] # initialise new current position
				current_state = binpart[i] # initialise new current state
		if i == len(binpart)-1: # conclude and append last partition
			current_part[1] = i+1
			if current_state == 0:  # append concluded partition to partitions_intron or partitions_exon, depending on current_state
				partitions_intron.append(current_part)
			else:
				partitions_exon.append(current_part)
	

	# write RAxML style partition file
	with open(gene+"_part.txt", "w") as partfile:
		print("DNA, intron = " + ", ".join(["-".join([str(j) for j in i]) for i in partitions_intron]), file=partfile)
		print("DNA, exon = " + ", ".join(["-".join([str(j) for j in i]) for i in partitions_exon]), file=partfile)
# write alignment without exon sequences
with open(gene+"_clean.fasta", "w") as al:	
	SeqIO.write(sequences, al, "fasta")