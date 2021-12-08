import numpy as np
from Bio import SeqIO
import fasta_to_seq as fts
import count_from_seqs as cfs

#class for mutating storing mutated and unmutated data on single organism
#input: single FNA file
#output: mutated FNA file and metadata
class mut_org_from_fasta():
	def __init__(self, fna_file, kmer_to_idx, r1=0):
		self.fna_file = fna_file
		self.r1 = r1
		self.k = len(list(kmer_to_idx.keys())[0])
		self.raw_seqs = fts.fasta_to_str_seq(self.fna_file)
		if r1 > 0:
			self.mutated_seqs = self.mutate_fasta_seqs()
		else:
			self.mutated_seqs = self.raw_seqs
		self.mut_kmer_ct, self.n_kmers_in_seqs = cfs.count_from_seqs(self.k,kmer_to_idx,self.mutated_seqs)

	def mutate_fasta_seqs(self):
		mut_seqs = []
		for seq in self.raw_seqs:
			seq = self.get_mutated_seq(seq)
			mut_seqs.append(seq)
		return mut_seqs

	def get_mutated_seq(self, seq):
		L = len(seq)
		mut_seq_list = list(seq)
		mut_flag = np.random.binomial(1,self.r1,L)
		mut_indices = np.nonzero(mut_flag)
		mut_indices = mut_indices[0]
		nucleotides = {'A','C','G','T'}
		for mut_idx in mut_indices:
			possible_muts = list(nucleotides.difference({seq[mut_idx]}))
			mut_seq_list[mut_idx] = np.random.choice(possible_muts)
		mut_seq = ""
		mut_seq = mut_seq.join(mut_seq_list)
		return mut_seq