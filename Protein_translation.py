#!usr/bin/env py

print ('starting installing packages \n \n ')

# Installing required packages
import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from collections import OrderedDict
import pandas as pd 
import re

# Table =1 means taking the human amino acid table
table = 1
min_pro_len = 10

# Create a counter that will help print only the phase1 header in the fasta file
count = 0 

# Open a file that will get all the samples with different length proteins
g = open("/proj/snic2020-2-10/private/Analyses/imke/proteomics/diversity_in_proteins/5_protein_fasta/differing_samples.fasta", "a")

# Open text file with data about genes
my_file = pd.read_csv("/proj/snic2020-2-10/private/Analyses/imke/proteomics/diversity_in_proteins/scripts/genes_set3.txt",sep='\t', header=None, index_col=0)
same_file = open("/proj/snic2020-2-10/private/Analyses/imke/proteomics/diversity_in_proteins/scripts/genes_set3.txt")
for line in same_file.readlines():
	theline = line.strip().split("\t") 
	gene =theline[0]
	starting_10_aa = (my_file.loc[gene][3])

	# Looping over all sequences in fasta file
	for seq_record in SeqIO.parse("/proj/snic2020-2-10/private/Analyses/imke/proteomics/diversity_in_proteins/4_dna_fasta/{}.fasta".format(gene), "fasta"):
		sequence = seq_record.seq
		simple_list=[]	
		sequence_id = seq_record.id
		individual_to_be = sequence_id.split("_")
		print(individual_to_be)
		individual = individual_to_be[0]+"_"+individual_to_be[1]
		print(individual)

		# Creating all possible proteins
		for strand, nuc in [(+1, sequence), (-1, sequence.reverse_complement())]:
			for frame in range(3):
				length = 3 * ((len(sequence)-frame) // 3)
				for pro in nuc[frame:frame+length].translate(table).split("*"):
					simple_list.append([str(pro),len(pro)])
					df=pd.DataFrame(simple_list,columns=['Sequence','Length'])

		# Sort all the proteins based on their length
		sorted_df = df.sort_values("Length", ascending = False)
		pd.set_option('display.max_colwidth', None)
		longest_protein_df = sorted_df['Sequence'].head(1)
		longest_protein = longest_protein_df.values[0]
		def_longest_protein = re.search(starting_10_aa+'(\w*)',longest_protein)
		if starting_10_aa in longest_protein:
			x = def_longest_protein.group()
			print(x)
			if len(x) != (my_file.loc[gene][2]):
				print("The length of the protein is not as in the genes_set3.txt file, it is {}".format(len(x)))
				g.write(str(gene+"\t"+individual+"\n"))
		else:
			x = ((my_file.loc[gene][2]) * "-")
			print(x)
		# Writing the longest protein to file
		f = open("/proj/snic2020-2-10/private/Analyses/imke/proteomics/diversity_in_proteins/5_protein_fasta/{}_protein.fasta".format(gene), "a")
                
		# Create a counter that will help print only the phase1 header in the fasta file
		count +=1
		if (count % 2) ==0:
			f.write(str(">"+individual+"_phase2"+"\n"+x+"\n"))
		else:
			f.write(str(">"+individual+"_phase1"+"\n"+x+"\n"))
		f.close()
same_file.close()
