import sys
import os
import pickle

species_name = sys.argv[1]
sample_name = sys.argv[2]
species_name = species_name.strip('.snps.gz')

#species_name = "Prevotella_copri"
#sample_name = "SRR7658580"
Nmax = 4
Lmax = 20000


output_directory = os.path.expanduser("/u/project/ngarud/wrshoema/snv_prevalence/data/strainfinder_output/")


#####FIRST CHECK IF SPECIES/SAMPLE HAS ALREADY BEEN PROCESSED#####
filename="/u/project/ngarud/wrshoema/snv_prevalence/data/strainfinder_output/strain_dict.pickle"
with open(filename, 'rb') as handle:
    strain_dict = pickle.load(handle)
if species_name in strain_dict.keys():
    if sample_name in strain_dict[species_name]["Africa"]:
        sys.stdout.write("passed")
        sys.exit()
    if sample_name in strain_dict[species_name]["North America"]:
        sys.stdout.write("passed")
        sys.exit()





####################
## You'll need to modify this script so the input/output not only contains what species to run strainfinder on, but also what host/sample to run it on.
###################
if not os.path.isdir( "%s%s/" % (output_directory, species_name)):
    os.makedirs("%s%s/" % (output_directory,species_name))

filename_prefix = "%s%s/%s" % (output_directory, species_name, sample_name)

#os.system('python create_StrainFinderInput.py -o %s --species %s --sample %s -L %d' % (output_directory, species_name, sample_name, Lmax))
#sys.exit(0)
for N in xrange(1,Nmax+1):

    #alignment_filename = filename_prefix+".strainfinder.p"
    alignment_filename = "/u/project/ngarud/wrshoema/snv_prevalence/data/strainfinder_input/%s/%s.pkl" % (species_name, sample_name)
    em_filename = filename_prefix+(".%d.em.cpickle" % N)
    log_filename = filename_prefix+(".%d.log.txt" % N)
    otu_filename = filename_prefix+(".%d.otu.txt" % N)

    os.system('python StrainFinder.py --aln %s -N %d --max_reps 10 --dtol 1 --ntol 2 --max_time 20000 --converge --em %s --em_out %s --otu_out %s --log %s --n_keep 2 --force_update --merge_out --msg' % (alignment_filename,N,em_filename, em_filename,otu_filename,log_filename))
