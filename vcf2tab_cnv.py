###############################################
# NAME : vcf2tab_cnv.py
# Date: 11/10/2023
version = "1.0"
###############################################


import argparse
import math
import sys
import pandas as pd
import os
from loguru import logger

def is_positive(number, SAMPLE):
	"""
	Questa funzione prende in input un numero 
	e restituisce True se il numero è positivo, False altrimenti.
	"""
	if number >= 0:
		return True
	else:
		n = open('negative_FC.log', 'a')
		n.write('[WARINIG] the sample '+SAMPLE+' has a Fold change in CNV with negative value\n')
		n.close()
		return False


def vcf_to_table(vcf_file, table_file, SAMPLE, MODE):

	if os.path.exists(table_file):
		MODE="a"
	else:
		MODE="w"
	with open(vcf_file, 'r') as vcf, open(table_file, MODE) as table:
		if MODE == 'a':
			pass
		else:
			table.write('ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n')
		SAMPLE=SAMPLE.split(".")[0] #remove .bam

		for line in vcf:

			if line.startswith("##fileformat"):
				version = (line.split("=")[-1]).strip()
				#logger.info(f"Analyzing {version} version")

			# Skip commented lines
			if line.startswith('##'):
				continue

			if line.startswith("#"):
				fields_names = line.split("\t")
				fields_names = list(map(lambda x: x.strip(), fields_names))
				chrom_position = fields_names.index("#CHROM")
				info_position = fields_names.index("INFO")
				start_position = fields_names.index("POS")
				qual_position = fields_names.index("QUAL")
				format_infos_position = fields_names.index("FORMAT")
				format_position = fields_names.index(SAMPLE)
				
				continue

			# Split the line by tabs
			fields = line.strip().split('\t')
			# Extract the data we want to keep
			chrom = fields[chrom_position].strip('chr')
			start = fields[start_position]
			infos = fields[info_position].split(";")
			end = [info.split("=")[-1] for info in infos if "END" in info][0]
			qual = fields[qual_position]
			info = fields[format_position]
		

			# Id = fields[2]
			# position=Id.split(":")[-1]
			# start=position.split("-")[0]
			# end=position.split("-")[1]

			
			# info = fields[7].split(';')
			# if len(info) == 2:
			# # 	end = info[0].split('=')[1]
			#  	gene =info[1].split('=')[1]
			# else:
			# # 	end = info[1].split('=')[1]
			# 	gene =info[2].split('=')[1]
				
			#fc = float(fields[9])

############### CASISTICHE POSSIBILI DA RICONTROLLARE CON LUCIANO ################################
			format = fields[format_infos_position]
			if version == "VCFv4.1":
				if format == "FC":
					fc = fields[-1]
			elif version == "VCFv4.2":
				format_infos = format.split(":")
				fc_position = format_infos.index("SM")
				sample_info = info.split(":")
				fc = sample_info[fc_position]
			else:
				version_number = version.split("v")[-1]
				logger.critical(f"The VCF file for sample {SAMPLE} is in an unsupported version (v{version_number}). Supported versions: [4.1, 4.2]")
				sys.exit()

			if fc != ".":
				fc = float(fc)
				if is_positive(fc, SAMPLE):
					log2fc = math.log(fc,2)
				else:
					fc = 0.0001
					log2fc = math.log(fc,2)
			else:
				continue
################ CHIEDERE A LUCIANO COSA FARE SE FC = . ############################

			# check negative Fold change values
			# if a negative fold chenge is found
			# the Fold change value is changed in 0.0001
			# if is_positive(fc, SAMPLE):
			# 	log2fc = math.log(fc,2)
			# else:
			# 	fc = 0.0001
			# 	log2fc = math.log(fc,2)
			# Write the data to the table file
			# segmentated data example
			# ID<TAB>chrom<TAB>loc.start<TAB>loc.end<TAB>num.mark<TAB>seg.mean
			table.write(f'{SAMPLE}\t{chrom}\t{start}\t{end}\t{qual}\t{log2fc}\n')
			#print(f'{SAMPLE}\t{chrom}\t{start}\t{end}\t{qual}\t{fc}\n')


def vcf_to_table_fc(vcf_file, table_file, SAMPLE, MODE):
	if os.path.exists(table_file):
		MODE="a"
	else:
		MODE="w"
	with open(vcf_file, 'r') as vcf, open(table_file, MODE) as table:
		if MODE == 'a':
			pass
		else:
			table.write('ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\tgene\tdiscrete\n')
		SAMPLE=SAMPLE.split(".")[0]
		
		for line in vcf:
			if line.startswith("##fileformat"):
				version = (line.split("=")[-1]).strip()
				#logger.info(f"Analyzing {version} version")

			if line.startswith('##'):
				continue

			if line.startswith("#"):
				fields_names = line.split("\t")
				fields_names = list(map(lambda x: x.strip(), fields_names))
				chrom_position = fields_names.index("#CHROM")
				info_position = fields_names.index("INFO")
				start_position = fields_names.index("POS")
				qual_position = fields_names.index("QUAL")
				format_position = fields_names.index("FORMAT")
				alt_position = fields_names.index("ALT")
				sample_position = fields_names.index(SAMPLE)
				
				continue

			fields = line.strip().split('\t')

			chrom = fields[chrom_position].strip('chr')
			start = fields[start_position]
			infos = fields[info_position].split(";")
			end = [info.split("=")[-1] for info in infos if "END" in info][0]
			qual = fields[qual_position]
			alt = fields[alt_position]
			format = fields[format_position]
			sample_infos = fields[sample_position]

			if version == "VCFv4.1":
				gene = [info.split("=")[-1] for info in infos if "ANT" in info][0]
				end = [info.split("=")[-1] for info in infos if "END" in info][0]
				if format == "FC":
					fc = fields[-1]
			elif version == "VCFv4.2":
				gene = [info.split("=")[-1] for info in infos if "SEGID" in info][0]
				format_infos = format.split(":")
				sample_infos = sample_infos.split(":")
				fc_position = format_infos.index("SM")
				fc = sample_infos[fc_position]
			

			if fc != ".":
				fc = float(fc)
				if not is_positive(fc, SAMPLE):
					fc = 0.0001
			else:
				continue
			
			########################
			# manage discrete data #
			########################

			if alt == "<DUP>":
				discr = "2"
			elif alt == "<DEL>":
				discr = "-2"
			else:
				discr = "0"
			
			# Write the data to the table file
			table.write(f'{SAMPLE}\t{chrom}\t{start}\t{end}\t{qual}\t{fc}\t{gene}\t{discr}\n')




def load_table(file_path):
	"""
	Load a table from a file into a Pandas DataFrame object.

	Parameters:
	file_path (str): path to the file containing the table.

	Returns:
	pandas.DataFrame: the loaded table.
	"""
	df = pd.read_csv(file_path, sep='\t',header=0)
	return df


def main(INPUT, OUTPUT, SAMPLE, MODE):

	vcf_to_table(INPUT, OUTPUT, SAMPLE, MODE)
	vcf_to_table_fc(INPUT, OUTPUT, SAMPLE, MODE)

if __name__ == '__main__':

	# parse arguments
	parser = argparse.ArgumentParser(description="Get a vcf file, generate a table",
											epilog="Version: 1.0\n\
											Author: Luciano Giaco'\n\
											email: luciano.giaco@policlinicogemelli.it")

    # arguments
	parser.add_argument('-i', '--input', help="<input.vcf>\
											VCF file for CNV",
											required=True)
	parser.add_argument('-o', '--output', help="<output-file.tab>\
											file path of the Table output",
											required=True)

	parser.add_argument('-s', '--sample', help="<sample name>\
											Name of the sample",
											required=True)
	parser.add_argument('-m', '--mode', help="Select the writing mode: write or append",
											choices=['w', 'a'],
											default='w')


	args = parser.parse_args()
	INPUT = args.input
	OUTPUT = args.output
	SAMPLE = args.sample
	MODE = args.mode

	main(INPUT, OUTPUT, SAMPLE, MODE)
