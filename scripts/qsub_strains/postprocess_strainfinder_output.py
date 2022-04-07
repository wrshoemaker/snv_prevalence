from StrainFinder import *
import cPickle
import sys
import numpy
import parse_midas_data
import parse_timecourse_data
from math import fabs
import os
import pickle

##############################
## this script will need to be modified so that it is also aware of the host/sample
##############################

species_name = sys.argv[1]
sample_name = sys.argv[2]
species_name = species_name.strip('.snps.gz')


#species_name = "Prevotella_copri_61740"
#sample_name = "SRR7658580"
Nmax = 4
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

sys.stderr.write("species name is...%s \n" % species_name)



filename_prefix = "%s%s/%s" % (output_directory, species_name, sample_name)

# Get best (min) AIC version
# Get filenames of EM objects
fns = [filename_prefix+(".%d.em.cpickle" % N) for N in range(1,Nmax+1)]

#species_name+sample_name+".strainfinder.locations.p"

# Load location object
#if not os.path.exists(filename_prefix+".strainfinder.locations.p"):
#    f = open(output_directory+'%s/final_strain_freq.txt' % species_name, 'a')
#    f.write(str( "Location file for sample %s is missing" % sample_name + '\n'))
#    f.close()

#snp_locations = cPickle.load(open(filename_prefix+".strainfinder.locations.p",'rb'))

# Load EM objects
ems = [cPickle.load(open(fn, 'rb')) for fn in fns]

# Get the best AIC in each EM object
bics = numpy.array([em.select_best_estimates(1)[0].aic for em in ems])
print "all BIC scores: ", bics
print "BIC scores:", bics.min()-bics

# Select EM with the minimum AIC
em = ems[numpy.argmin(bics)]
#em = ems[1]

# M = timepoint, L = along genome, 4 = alleles
input_alignment = em.data.x # alignment data, dim = (M x L x 4)

best_estimate = em.select_best_estimates(1)[0]

strain_genotypes = best_estimate.p # genotypes of first estimate, dim = (N x L x 4)
strain_freqs = best_estimate.z # (M x N)

#output_alignment = numpy.dot(strain_freqs, strain_genotypes)

strain_genotypes_derived = strain_genotypes[:,:,0]
strain_genotypes_ancestral = strain_genotypes[:,:,1]

strain_counts = strain_genotypes_derived.sum(axis=0)
strain_counts = numpy.fmin(strain_counts,strain_genotypes.shape[0]-strain_counts)
#print "SFS strains"

singletons = (strain_counts==1)

#for k in xrange(0,strain_genotypes.shape[0]):
    #print k, (strain_counts==k).sum()

output_alignment = numpy.einsum('ij,jkl', strain_freqs, strain_genotypes)

#print input_alignment.shape
#print output_alignment.shape

print "Best # of strains:", em.select_best_estimates(1)[0].z.shape[1]
print em.select_best_estimates(1)[0].z # frequencies of first estimate, dim = (M x N)

f = open(output_directory+'%s/final_strain_freq.txt' % species_name, 'a')
#f.write(str("all BIC scores: ") + str(bics))
f.write(str( "Best # of strains for sample %s: " % sample_name + str(em.select_best_estimates(1)[0].z.shape[1])+ ' with frequencies: '+ str(em.select_best_estimates(1)[0].z) + '\n'))
f.close()
