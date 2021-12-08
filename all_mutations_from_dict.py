import os
import numpy as np
from Bio import SeqIO
import glob
import warnings
import MinHash as mh
import itertools
from pathlib import Path
import mut_org_from_fasta as moff
import re

#inputs: org_dict object, sample abundance (list), sample mutation rates (list)
#output: all_mutations_from_dict object, contains multiple mutated samples and aggregated kmer frequency vector
class all_mutations_from_dict():
	def __init__(self, org_dict, rel_abundance,r1,total_kmers = pow(10,9), rnd = True):
		self.OD = org_dict
		self.mut_kmer_ct = np.zeros((len(list(self.OD.kmer_to_idx.keys())),))
		self.rel_abundance = rel_abundance
		self.r1 = r1
		self.total_kmers = total_kmers
		self.mut_orgs = []

		for file_idx, fasta_file in enumerate(self.OD.fasta_files[0:len(r1)]):
		    curr_mut_org = moff.mut_org_from_fasta(fasta_file, self.OD.kmer_to_idx, self.r1[file_idx])
		    self.mut_orgs.append(curr_mut_org)
		    self.mut_kmer_ct = self.mut_kmer_ct + self.total_kmers*self.rel_abundance[file_idx]*curr_mut_org.mut_kmer_ct
		if rnd:
			self.mut_kmer_ct = np.round(self.mut_kmer_ct)