
0. Use write_genome_file.py to quickly create a text file with names of all FNA files from a directory for easy use with MinHash.

1. Run MinHash to create a count estimator database (.h5) file.

2. Create the dictionary (A matrix) using org_dict.py, which takes as input the path to the h5 file in the previous step.

3. Create simulated mutated data using all_mutations_from_fasta.py: create a all_mutations_from_fasta object which has 4 main inputs.
	-org_dict: dictionary object from previous step
	-rel_abundance: list of abundance of the organisms in the dictionary in the sample, ex [0.1, 0.4, 0.5]
	-r1: list of mutation rates for simulated data 

4. Create a frequency_estimator object using freq_est_from_mut.py. Takes as input an all_mutations_from_dict object. Runs the linear program through lp_solver.py and computes the estimated frequency as frequency_estimator.freq_est

5. mut_test_script.py runs through these steps and prints output of step 4.