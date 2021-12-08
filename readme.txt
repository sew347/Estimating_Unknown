
0. Use write_genome_file.py to quickly create a text file with names of all FNA files from a directory for easy use with MinHash.

1. Run MinHash to create a count estimator database (.h5) file.

2. Create the dictionary (A matrix) using org_dict.py, which takes as input the path to the h5 file in the previous step.

3. Create simulated mutated data using all_mutations_from_fasta.py: create a all_mutations_from_fasta object which has 4 main inputs.
	-org_dict: dictionary object from previous step
	-rel_abundance: list of abundance of the organisms in the dictionary in the sample, ex [0.1, 0.4, 0.5]
	-r1: list of mutation rates for simulated data 

4. Create a frequency_estimator object using freq_est_from_mut.py. Takes as input an all_mutations_from_dict object. Runs the linear program through lp_solver.py and computes the estimated frequency as frequency_estimator.freq_est

5. mut_test_script.py runs through these steps and prints output of step 4.

Example command and expected output below. In this example, the first two organisms are recovered accurately, while the last is rejected due to having a mutation rate above the threshold of 0.05.

stephenwhite$ python mut_test_script.py -a 0.3 0.2 0.5 -r1 0.01 0.02 0.07 -seed 10

Academic license - for non-commercial use only - expires 2022-08-11
Using license file /Users/stephenwhite/gurobi.lic
Recovered frequencies at 0.05 mutation threshold:
[0.3005 0.2002 0.    ]
