import os
import argparse

def write_genome_file(filepath, genome_dir, N = None):
	with open(filepath,'w+') as fid:
		org_names = os.listdir(genome_dir)
		if N is not None:
			org_names = org_names[:N]
		for idx, file_name in enumerate(org_names):
			if file_name.endswith('.fna'):
			    fid.write(os.path.abspath(genome_dir+'/'+file_name))
			    if idx < len(org_names)-1:
			    	fid.write('\n')
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('genome_dir')
	parser.add_argument('out_file')
	args = parser.parse_args()
	genome_dirname = os.path.abspath(args.genome_dir)
	if not os.path.isdir(genome_dirname):
		raise Exception("Genome dir %s does not exist." % genome_dirname)
	out_file = os.path.abspath(args.out_file)
	write_genome_file(out_file, genome_dirname)

if __name__ == "__main__":
	main()
