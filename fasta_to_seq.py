from Bio import SeqIO
import re

def fasta_to_str_seq(fasta):
	notACTG = re.compile('[^ACTG]')
	fasta_seqs = SeqIO.parse(open(fasta),'fasta')
	seqs = []
	for seq_iter in fasta_seqs:
		seq_str = str(seq_iter.seq).upper()
		seq_split_ACTG = notACTG.split(seq_str)
		for seq in seq_split_ACTG:
			seqs.append(seq)
	return seqs