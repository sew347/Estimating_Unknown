import os
import numpy as np
from Bio import SeqIO
import glob
import warnings
import MinHash as mh
import itertools
from pathlib import Path
import pickle
import re
import fasta_to_seq as fts
import count_from_seqs as cfs

#takes in count estimator file and creates dictionary matrix and associated metadata
#input: 
class org_dict():
	def __init__(self, db_file, N = None, filename = None, filepath = None):
		self.dict_files_from_db(db_file, N=N)
		if filename is not None:
			self.filename = filename
		else:
			self.filename = self.dict_filename(filepath,N=N)
		if self.filename is not None:
			with open(self.filename, 'wb') as f:
				pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

	def dict_files_from_db(self,db_file, N=None):
		metadata = mh.get_info_from_single_hdf5(db_file)
		if N is None:
			N = len(metadata.file_names)
		genome_files = (metadata.file_names)[0:N]
		genome_data = mh.import_multiple_from_single_hdf5(db_file, genome_files)
		temp = self.kmer_union_from_count_estimators(genome_data)
		idx_to_kmer = temp[0]
		kmer_to_idx = temp[1]
		self.k = len(idx_to_kmer[0])
		self.n_kmers = len(genome_data[0]._kmers)
		self.fasta_files = genome_files
		self.idx_to_kmer = idx_to_kmer
		self.kmer_to_idx = kmer_to_idx
		dict_matrix = self.matrix_from_fasta_files()
		self.dictionary = dict_matrix

	def matrix_from_fasta_files(self):
		n_kmers = len(self.kmer_to_idx.keys())
		N = len(self.fasta_files)
		org_dict = np.zeros((n_kmers,N))
		for curr_idx, curr_fasta in enumerate(self.fasta_files):
			curr_seqs = fts.fasta_to_str_seq(curr_fasta)
			org_dict[:,curr_idx], temp = (cfs.count_from_seqs(self.k, self.kmer_to_idx, curr_seqs))
		return org_dict

	def kmer_union_from_count_estimators(self,count_estimators):
		idx_to_kmer = []
		kmer_to_idx = {}
		org_kmers = []
		for ce in count_estimators:
			org_kmers.append(ce._kmers)
			idx_to_kmer.extend(ce._kmers)
		idx_to_kmer = list(set(idx_to_kmer))
		for idx, kmer in enumerate(idx_to_kmer):
			kmer_to_idx[kmer] = idx
		return [idx_to_kmer, kmer_to_idx]

	def dict_filename(self,filepath, N=None):
		if os.path.isfile(filepath):
			dict_dir = os.path.dirname(filepath)
		elif os.path.isdir(filepath):
			dict_dir = filepath
		else:
			raise ValueError('Argument is not a file or directory.')
		if N is not None:
			add = '_N'+str(N)
		else:
			add = ''
		dict_files = dict_dir+'/org_dict'+add+'.pkl'
		return dict_files
