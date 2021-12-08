import os
import numpy as np
from Bio import SeqIO
import MinHash as mh
import re

def get_kmer_count(k, kmer_to_idx, seq):
	L = len(seq);
	n_kmers = len(kmer_to_idx.keys())
	count = np.zeros((n_kmers,))
	for i in range(L-k+1):
		kmer = seq[i:(i+k)]
		if kmer in kmer_to_idx:
			count[kmer_to_idx[kmer]] += 1
	return count

def count_from_seqs(k, kmer_to_idx, seqs):
	keys = kmer_to_idx.keys()
	n_kmers = len(keys)
	n_kmers_in_seqs = 0
	count = np.zeros((n_kmers,))
	for seq in seqs:
		n_kmers_in_seqs += (len(seq) - k + 1)
		count = count + get_kmer_count(k,kmer_to_idx,seq)
	if np.sum(count) == 0:
		warnings.warn('Count vector is empty.');
	count = count/n_kmers_in_seqs
	return count, n_kmers_in_seqs