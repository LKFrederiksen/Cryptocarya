#!/usr/bin/python3

'''
------------------------------------------------------------------------------------------------------------------------
This workflow is used in a Workflow to estimate coverage.
------------------------------------------------------------------------------------------------------------------------
This code is a variant of coverage.py by Wolf Eiserhardt
Edited by Sarah E.K. Kessel
------------------------------------------------------------------------------------------------------------------------
'''

import os, argparse, subprocess
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("sample")
parser.add_argument("directory_in")
parser.add_argument("directory_out")
parser.add_argument("directory_wrk")
args = parser.parse_args()
sample = str(args.sample)
directory_in = str(args.directory_in)
directory_out = str(args.directory_out)
directory_wrk = str(args.directory_wrk)

# depth required to KEEP (i.e. anything <trshld will be discarded)
trshld = 2

# Go to working directory
cmd = 'cd '+directory_wrk
subprocess.call(cmd,shell=True)

# Get all subdirectories in the current working directory. these are the loci recovered by hybpiper
loci = next(os.walk(sample))[1]
sequences = {}

##### Try by using 'a' instead of 'w' #####
# 90% of worst supercontig partitions removed.
loci_supercontig = ["4951", "4992", "5116", "5131", "5162", "5168", "5328", "5335", "5427", "5554", "5594", "5614", "5620", "5770", "6016", "6036", "6064", "6098", "6139", "6318", "6363", "6488", "6500", "6563" "6913", "7029", "7111"]
for supercontig in loci_supercontig:
        pth = sample+'/'+supercontig+'/'+sample+'/sequences/intron/'+supercontig+'_supercontig.fasta'

        if os.path.isfile(pth):
                for record in SeqIO.parse(pth, "fasta"):
                        record.id = record.id+'_'+supercontig
                        #sequences.append(record)
                        sequences[record.id] = record

with open(directory_out+sample+'.fasta', "w") as outfile:
        SeqIO.write(list(sequences.values()), outfile, "fasta")

print(sample+'.fasta generated')

# exon
loci_exon = ["6460", "6384", "7313", "5858", "6782", "6420", "5032", "5318", "6883", "5333", "7572", "6494", "5639", "5339", "4954", "5899", "6909", "4802", "5716", "6373", "5404", "6462", "6570", "5271", "5699", "6528", "5849", "6639", "6992", "5463", "5426", "5893", "5428", "6914", "4527", "6961", "6544", "5449", "5034", "6051", "7136", "5702", "7577", "5454", "6533", "6457", "5163", "5894", "6924", "4932", "5348", "6447", "7371", "5421", "5280", "5870", "6496", "5791", "5945", "5634", "7363", "6072", "5038", "5343", "5913", "5596", "7324", "6458", "5857", "5816", "5942", "5981", "6295", "6641", "4942", "5815", "5406", "6689", "5326", "4757", "6779", "4793", "6954", "7024", "6652", "6848", "6933", "6003", "6978", "7279", "5926", "5502", "6068", "6404", "5477", "6738", "6631", "5138", "6649", "5977", "6432", "6450", "6459", "6216", "6968", "7067", "7021", "6376", "5177", "7141", "7336", "6865", "7174", "6732", "4691", "4724", "5910", "6947", "5562", "4848", "6175", "6282", "6038", "6412", "6483", "6667", "5464", "5264", "5664", "6860", "6882", "5865", "5104", "5670", "5220", "6733", "6636", "5703", "5551", "5950", "6389", "6559", "5200", "7241", "5644", "6572", "5018", "5721", "7331", "6454", "5366", "6685", "6439", "6825", "5257", "5528", "7367", "5296", "5206", "4890", "6000", "6298", "6492", "5918", "5843", "5273", "5802", "7325", "7333", "6299", "6378", "5536", "5513", "4889", "6527", "5933", "6004", "5968", "6550", "6284", "6383", "5531", "6148", "6226", "5299", "6274", "6538", "5744", "5489", "6176", "6797", "6029", "5469", "5188", "5398", "6601", "6227", "5949", "5578", "5980", "5990", "7194", "6398", "6050", "5841", "5840", "6532", "6198", "6258", "6303", "5821", "6552", "6265", "5960", "6238", "6320", "6875", "5919", "6717", "5090", "6792", "4796", "7628", "6979", "6746", "6407", "5944", "6854", "5460", "4806", "6366", "6405", "6507", "6048", "5822", "6946", "5866", "6859", "5123", "6026", "5974", "6962", "5304", "6164", "5772", "6056", "6785", "5853", "6114", "6128", "7128", "6034", "4471", "6130", "6660", "6958", "4744", "7135", "6506", "6620", "6221", "6713", "6393", "5354", "5660", "7273", "6526"]
for exon in loci_exon:
        pth = sample+'/'+exon+'/'+sample+'/sequences/FNA/'+exon+'.FNA'
                
        if os.path.isfile(pth):
                for record in SeqIO.parse(pth, "fasta"):
                        record.id = record.id+'_'+exon
                        #sequences.append(record)
                        sequences[record.id] = record
           
with open(directory_out+sample+'.fasta', "a") as outfile:
        SeqIO.write(list(sequences.values()), outfile, "fasta")

print(sample+'.fasta with exons generated')


# BWA index targets
cmd = 'bwa index '+directory_out+sample+'.fasta'
subprocess.call(cmd,shell=True)
print(sample+'.fasta indexed')

# BWA mem paired reads and @HD tag
cmd = 'bwa mem '+directory_out+sample+'.fasta '+directory_in+sample+'_1P.fastq '+directory_in+sample+'_2P.fastq | samtools view -b -h -o '+directory_out+sample+'.bam'
subprocess.call(cmd,shell=True)
print('paired reads mapped to '+sample+'.fasta')

# BWA mem unpaired reads
cmd = 'bwa mem '+directory_out+sample+'.fasta '+directory_in+sample+'_UN.fastq | samtools view -b -h -o '+directory_out+sample+'_up.bam'
subprocess.call(cmd,shell=True)
print('unpaired reads mapped to '+sample+'.fasta')

# @HD
# cmd = 'bam polishbam --in '+directory_out+sample+'_no_up.bam --out '+directory_out+sample+'_up.bam --HD "@HD    VN:1.3 SO:coordinate"'
# subprocess.call(cmd,shell=True)
# cmd = 'bam polishbam --in '+directory_out+sample+'_no.bam --out '+directory_out+sample+'.bam --HD "@HD  VN:1.3 SO:coordinate"'
# subprocess.call(cmd,shell=True)
# print('@HD added')

# merge BAM files
cmd = 'samtools merge -f '+directory_out+sample+'_all.bam '+directory_out+sample+'.bam '+directory_out+sample+'_up.bam'
subprocess.call(cmd,shell=True)
print('BAMs merged')

# sort and index BAM files
cmd = 'samtools sort '+directory_out+sample+'_all.bam -o '+directory_out+sample+'_all_sorted.bam'
subprocess.call(cmd,shell=True)
cmd = 'samtools index '+directory_out+sample+'_all_sorted.bam'
subprocess.call(cmd,shell=True)
print('BAM indexed and sorted')

# remove duplicates
cmd = 'picard MarkDuplicates I='+directory_out+sample+'_all_sorted.bam O='+directory_out+sample+'_all_sorted_deduplicated.bam M='+directory_out+sample+'_marked_dup_metrics.txt REMOVE_DUPLICATES=true'
subprocess.call(cmd,shell=True)
print('reads deduplicated for sample '+sample)

# calculate coverage # Removed _deduplicated because the above did not work.
cmd = 'samtools depth '+directory_out+sample+'_all_sorted.bam > '+directory_out+sample+'.cov'
subprocess.call(cmd,shell=True)
print('coverage calculated for sample '+sample)

# define function to replace nth position of sequence with N
def n2N(sqnc, pstn):
        sqnc = list(sqnc)
        sqnc[int(pstn)-1] = "N"
        return "".join(sqnc)

# process coverage
with open(directory_out+sample+'.cov', "r") as covfile:
        for line in covfile:
                line = line.strip()
                LINE = line.split("\t")
                if int(LINE[2]) < trshld:
                        sequences[LINE[0]].seq = n2N(sequences[LINE[0]].seq,LINE[1])

# remove unnecessary leading and trailing Ns
for nm in sequences.keys():
        sequences[nm].seq = sequences[nm].seq.strip("N")
        if isinstance(sequences[nm].seq, str):
                sequences[nm].seq = Seq(sequences[nm].seq)

print('coverage trimming completed, keeping only positions with coverage of '+str(trshld)+' or above')

# write outfile
with open(directory_out+sample+'_trimmed.fasta', "w") as outfile:
        SeqIO.write(list(sequences.values()), outfile, "fasta")
print('trimmed seqs written to '+sample+'_trimmed.fasta')

# remove unnecessary files
#cmd = "find ../coverage -type f ! -name '*.fasta' -delete"
#subprocess.call(cmd, shell=True)

