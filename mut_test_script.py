import numpy as np
import org_dict as od
import all_mutations_from_dict as amfd
import lp_solver as lps
import freq_est_from_mut as ferm
import argparse

def main():
	parser = argparse.ArgumentParser(description="This script runs random tests on a given organism dictionary.", \
	 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('-a,','--abundance', help='Relative abundance of organisms in simulated sample', default=0.05, type = float, nargs = '*')
	parser.add_argument('-r1', help='Mutation rate of organisms in simulated sample', default=0.05, type = float, nargs = '*')
	parser.add_argument('-db_file', default = "files/databases/cdb_n1000_k31/cdb.h5", help='file with count estimators')
	parser.add_argument('-w', '--w', help='False negative discount weight', default=None, type = float)
	parser.add_argument('-seed', '--seed', help='Random seed', default=None, type = int)

	args = parser.parse_args()
	
	if args.seed is not None:
		np.random.seed(args.seed)

	db_file = args.db_file
	abundance = args.abundance
	r1 = args.r1
	N = len(abundance)
	if len(r1) != len(abundance):
		raise ValueError("r1 and abundance must be same length.")
	w = args.w

	OD = od.org_dict(db_file=db_file,N=N, filepath = "files/databases/cdb_n1000_k31")
	k = OD.k

	AM = amfd.all_mutations_from_dict(OD, abundance, r1)

	FE = ferm.frequency_estimator(AM,w=w)
	print("Recovered frequencies at " + str(0.05) + " mutation threshold:")
	print(np.round(FE.freq_est,4))

if __name__ == "__main__":
	main()