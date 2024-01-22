#!/usr/bin/env python3

from cyvcf2 import VCF
import numpy as np
import argparse
import re, sys, gc, os

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(
			prog="VCFpartitioner",
			description="""Multiallelic variants_tsv should be split and variants_tsv should be atomized (bcftools norm VCF.vcf.gz -m- -a -N)""",
			epilog="Written by Adrien Le Meur, v.0.1")

parser.add_argument('-i', type=str, nargs='+', required=True, help='a single sample VCF aligned on H37Rv genome')
parser.add_argument('-o', type=str, nargs='+', required=False, default = "./", help='output directory')

parser.add_argument('--informative', default = False, action='store_true', help='output only informative sites in MSA')
parser.add_argument('--bed', default = False, type=argparse.FileType('r'), required=False, help='bed file with column the annotation')

args = parser.parse_args()

informative = args.informative
input_vcf = args.i
input_bed = args.bed

output_directory = str(args.o[0])

if(args.o):
	output_directory = str(args.o[0])
	if not os.path.exists(output_directory):
		os.mkdir(output_directory)

all_samples = [ VCF(sample).samples[0] for sample in input_vcf ]


bed_dict = {}
if args.bed:
	print("Reading a bed file")
	for line in args.bed:
		line = line.rstrip('\n').split("\t")
		key = "_".join(line)
		bed_dict[key] = {"CHROM":line[0], "START":min(line[1],line[2]), "STOP":max(line[1],line[2]), "ANN":line[3]}

def search_annotation(POS, my_bed_dict, switch):
	if switch:
		annotation = []
		for i in my_bed_dict:
			if POS >= int(my_bed_dict[i]["START"]) and POS <= int(my_bed_dict[i]["STOP"]):
				annotation.append(my_bed_dict[i]["ANN"])
		if annotation == []:
			return("NONE")
		else:
			return("/".join(annotation))
	else:
		return("NONE")

#all_samples = all samples found in files

#SNP
#key = { POS_REF_ALT :	{POS, REF, { sample : VARIANT } }

#INDEL
#key = { POS_REF_ALT :	{POS, REF, ALT, { sample : 1|0 } }

# /!\ both dictionnaries are different /!\
#could be simplified for optimisation

MSA_dict = {}
MSA_dict['SNP'] = {}
MSA_dict['INDEL'] = {}

for sample in input_vcf:
	vcf = VCF(sample)

	sample_name = vcf.samples[0]

	for v in vcf:
		for variant in v.ALT:
			if v.FILTER == 'FAIL' : continue

			key = "_".join([str(v.POS), v.REF, v.ALT[0]])

			if v.is_snp or v.is_mnp:
				if key not in MSA_dict['SNP']:
					MSA_dict['SNP'][key] = {}
					MSA_dict['SNP'][key]['POS'] = v.POS
					MSA_dict['SNP'][key]['REF'] = v.REF

					for i in all_samples:
						if i == sample_name:
							#new entry in the dictonnary
							#until further notice, all other samples have the reference phenotype
							MSA_dict['SNP'][key][i] = variant
						else:
							MSA_dict['SNP'][key][i] = v.REF
				else:
					MSA_dict['SNP'][key][sample_name] = variant

			elif v.is_indel:
				if key not in MSA_dict['INDEL']:
					MSA_dict['INDEL'][key] = {}
					MSA_dict['INDEL'][key]['POS'] = v.POS
					MSA_dict['INDEL'][key]['REF'] = v.REF
					MSA_dict['INDEL'][key]['ALT'] = v.ALT[0]

					for i in all_samples:
						if i == sample_name:
							MSA_dict['INDEL'][key][i] = 1
						else:
							MSA_dict['INDEL'][key][i] = 0
				else:
					MSA_dict['INDEL'][key][sample_name] = 1

#second dictionnaries parse to check whether variants are informative
for i in MSA_dict['SNP']:
	pattern = {}
	for j in MSA_dict['SNP'][i]:
		if j != 'POS' and j != 'REF':
			if MSA_dict['SNP'][i][j] not in pattern:
				pattern[MSA_dict['SNP'][i][j]] = 1
			else:
				pattern[MSA_dict['SNP'][i][j]] += 1
	if(len(pattern.keys()) > 1 and min(pattern.values()) >= 2):
		MSA_dict['SNP'][i]['INFORMATIVE'] = 'TRUE'
	else:
		MSA_dict['SNP'][i]['INFORMATIVE'] = 'FALSE'

for i in MSA_dict['INDEL']:
	pattern = {}
	for j in MSA_dict['INDEL'][i]:
		if j != 'POS' and j != 'REF' and j != 'ALT':
			if MSA_dict['INDEL'][i][j] not in pattern:
				pattern[MSA_dict['INDEL'][i][j]] = 1
			else:
				pattern[MSA_dict['INDEL'][i][j]] += 1
	if(len(pattern.keys()) > 1 and min(pattern.values()) >= 2):
		MSA_dict['INDEL'][i]['INFORMATIVE'] = 'TRUE'
	else:
		MSA_dict['INDEL'][i]['INFORMATIVE'] = 'FALSE'

#MSA = { SAMPLE : { [BIN] : str() }, {[NUC] : str() }}
MSA = {}

for i in all_samples:
	MSA[i] = {}
	MSA[i]['NUC'] = ''
	MSA[i]['BIN'] = ''

if(args.o):
	variants_tsv = open(str(output_directory+"/variants.tsv"), "w")
else:
	variants_tsv = open("variants.tsv", "w")

previous_variant = {}

variants_tsv.write("START\tSTOP\tREF\tALT\tINFORMATIVE\tANN\n")

#parse both dictionnaries to print a list of variants
for i in MSA_dict['SNP']:
	variants_tsv.write(str(MSA_dict['SNP'][i]['POS'])+"\t")
	variants_tsv.write(str(MSA_dict['SNP'][i]['POS']+1)+"\t")
	variants_tsv.write(str(MSA_dict['SNP'][i]['REF'])+"\t")
	for j in MSA_dict['SNP'][i]:
		if j != 'POS' and j != 'REF' and j != 'INFORMATIVE':
			#append here
			if( informative is False or MSA_dict['SNP'][i]['INFORMATIVE'] == 'TRUE' ):
				MSA[j]['NUC']+=str(MSA_dict['SNP'][i][j])

			if(MSA_dict['SNP'][i][j] not in previous_variant and MSA_dict['SNP'][i][j] != MSA_dict['SNP'][i]['REF']):
				variants_tsv.write(str(MSA_dict['SNP'][i][j]))
				#get all polymorphisms
				previous_variant[MSA_dict['SNP'][i][j]]=1
	variants_tsv.write("\t"+str(MSA_dict['SNP'][i]['INFORMATIVE']))	
	variants_tsv.write("\t"+search_annotation(MSA_dict['SNP'][i]['POS'], bed_dict, input_bed))
	previous_variant = {}
	variants_tsv.write("\n")


for i in MSA_dict['INDEL']:
	variants_tsv.write(str(MSA_dict['INDEL'][i]['POS'])+"\t")
	variants_tsv.write(str(MSA_dict['INDEL'][i]['POS']+1)+"\t")
	variants_tsv.write(str(MSA_dict['INDEL'][i]['REF'])+"\t")
	variants_tsv.write(str(MSA_dict['INDEL'][i]['ALT'])+"\t")
	variants_tsv.write(str(MSA_dict['INDEL'][i]['INFORMATIVE'])+"\t")
	variants_tsv.write(search_annotation(MSA_dict['INDEL'][i]['POS'], bed_dict, input_bed)+"\n")
	for j in MSA_dict['INDEL'][i]:
		if j != 'POS' and j != 'REF' and j != 'ALT' and j != 'INFORMATIVE':
			#append here
			if( informative is False or MSA_dict['INDEL'][i]['INFORMATIVE'] == 'TRUE' ):
				MSA[j]['BIN']+=str(MSA_dict['INDEL'][i][j])

fasta = open(str(output_directory+"/nuc_MSA.fasta"), "w")
fastb = open(str(output_directory+"/binary_MSA.fasta"), "w")
fastm = open(str(output_directory+"/mixed_MSA.fasta"), "w")
partition = open(str(output_directory+"/raxml_partition.txt"), "w")

first_key = list(MSA.keys())[0]
partition.write("GTR+G, p1=1-"+str(len(MSA[first_key]['NUC']))+"\n"+"BIN, p2="+str(len(MSA[first_key]['NUC'])+1)+"-"+str(len(MSA[first_key]['NUC'])+len(MSA[first_key]['BIN']))+"\n")

for i in MSA:
	fastm.write(">"+i+"\n")#+"_"+'NUC_LENGTH/'+str(len(MSA[i]['NUC']))+'/BIN_LENGTH/'+str(len(MSA[i]['BIN']))+"/\n")
	fasta.write(">"+i+"\n")#+"_"+'NUC_LENGTH/'+str(len(MSA[i]['NUC']))+"/\n")
	fastb.write(">"+i+"\n")#+"_"+'BIN_LENGTH/'+str(len(MSA[i]['BIN']))+"/\n")

	fastm.write(MSA[i]['NUC']+MSA[i]['BIN']+"\n")
	fasta.write(MSA[i]['NUC']+"\n")
	fastb.write(MSA[i]['BIN']+"\n")

fastm.close()
fasta.close()
fastb.close()
