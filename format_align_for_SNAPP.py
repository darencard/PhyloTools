#!/usr/bin/env python

usage_line = """
format_align_for_SNAPP.py

Version 1.0 (17 September, 2014)
License: GNU GPLv2
To report bugs or erros, please contact Daren Card (dcard@uta.edu)
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

Script takes a sequence alignment and converts it to a binary nexus output that can be used in 
the program SNAPP. Input can be in any alignment format supported by the AlignIO BioPython package,
with phylip, nexus, and fasta being the most common formats encountered. User must specify the format
using the '--type' flag. User can also specify the output file name desired. The '--missing' flag
allows for missing data to be kept in the output. If not specified, any locus with missing data will be
excluded. The script also checks that loci are bi-allelic and ignores any locus that is not.

python format_align_for_SNAPP.py --input <input.ext. --output <output.nex> --type <nexus/fasta/phylip/etc>
[--missing]
"""

import os
import sys
import optparse
from Bio import AlignIO

usage = usage_line

parser = optparse.OptionParser(usage=usage)
parser.add_option("--input", action="store", type="string", dest="input", help="""The input alignment file""")
parser.add_option("--output", action="store", type="string", dest="output", help="""The file name for SNAPP nexus output [output.nex]""", default="output.nex")
parser.add_option("--type", action="store", type="string", dest= "type", help="""The type of input alignment file (nexus, phylip, fasta, etc.)""")
parser.add_option("--missing", action="store", type="string", dest="miss", help="""The amount of missing data allowed (0-1), after which a locus is excluded""")
#parser.add_option("--missing", action="store_true", dest="miss", help="""Specify flag to allow missing data""")

options, args = parser.parse_args()

def reformat_alignment(new_align):	
	alignment = AlignIO.read(options.input, options.type)
	for record in alignment:
		print record.id
#		print len(record)
		misscount = record.seq.count('N')
#		print misscount
		taxamiss = float(misscount)/len(record)
		print taxamiss
	new_alignment = list()
	for w in xrange(alignment.get_alignment_length()):
		bases = alignment[:, w] 	
		ncount = bases.count('N')
		length = float(len(bases))
		missing = ncount/length
		if missing < float(options.miss):
			uniqs = list()
			uniqs = list(set(bases))
			nuniqs = filter(lambda a: a != "N", uniqs) # Filter out N's to check for biallelic-ness
			if len(nuniqs) != 2: # Check for biallelic-ness
#				print "Skipping site {0} - not a biallelic SNP".format(w)
				nothing = 1
			else:
#				print "Including site {0}".format(w)
				allele1 = nuniqs[0]
				allele2 = nuniqs[1]
				new_snp = [0]*(len(bases))				
				for x in range(len(bases)):
					if bases[x] == allele1:
						new_snp[x] = 0
					elif bases[x] == allele2:
						new_snp[x] = 1
					elif bases[x] == "N":
						new_snp[x] = "?"			
				new_alignment.append(new_snp)
#		else:
#			print "Skipping site {0} - too much missing data".format(w)
#	print len(new_alignment)
	write_outfile(alignment, new_alignment)


def write_outfile(alignment, new_alignment):
	out = open(options.output, 'wb')
	out.write("#NEXUS\n\n")
	out.write("Begin data;\n")
	out.write("\tDimensions ntax={0} nchar={1};\n".format(len(alignment[:,1]), len(new_alignment)))
	out.write("\tFormat datatype=binary symbols=\"01\" gap=- missing=?;\n")
	out.write("\tMatrix\n")
	for x in range(len(alignment[:,1])): # Loop through individuals
		out.write(alignment[x].id+"\t")				
		for y in range(len(new_alignment)): # Loop over alignment columns
			out.write("{0}".format(new_alignment[y][x])) # Write nucleotides
			out.flush()
		out.write("\n")	
		out.flush()
	out.write("\t;\n")
	out.write("End;")
	out.close()
	

def main():
	new_align = list()
	reformat_alignment(new_align)
	
main()
	