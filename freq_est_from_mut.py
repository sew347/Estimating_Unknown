import all_mutations_from_dict as amfd
import lp_solver as lps
import org_dict as od
import numpy as np

#inputs: all_mutations_from_fasta object (simulated data)
#outputs: frequency estimator
class frequency_estimator():
	def __init__(self, all_mutations_from_dict, r1_thresh = 0.05, w = None):
		self.AM = all_mutations_from_dict
		self.k = self.AM.OD.k
		self.n_kmers = self.AM.OD.n_kmers
		self.r1_thresh = r1_thresh
		self.w = w if w is not None else self.estimate_w()
		self.LPS = lps.lp_solver(self.AM.OD.dictionary, self.AM.mut_kmer_ct, self.w, run_now = True)
		self.count_est = self.LPS.x_opt
		self.freq_est = self.LPS.x_opt/self.AM.total_kmers
		self.est_unk_pct = np.sum(self.freq_est)

	def estimate_w(self,est_n_orgs = 1000, p_val = 0.95, n_tests = 10000):
		prob = (1-self.r1_thresh)**self.k
		b = []
		for i in range(n_tests):
			b.append(min(np.random.binomial(self.n_kmers, prob, (est_n_orgs, 1))))
		min_est = np.quantile(b,p_val)
		print(min_est)
		w = min_est/(self.n_kmers-min_est)
		print(w)
		return w
