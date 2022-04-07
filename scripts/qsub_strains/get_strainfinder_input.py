import numpy as np
import pandas as pd
import sys
import os
import pickle

def get_counts(counts,depth_thresh):

    emp = []

    for i in range(counts.shape[0]):

        p = counts[i].split(',')
        elem = [int(e) for e in p]

        ## here, we look for sites that have more than 0 reads mapping
        ## to the second most common nucleotide e.g. polymorphic sites
        ## we append such sites to the list emp
        e_s = sorted(elem)
        s = sum(e_s)
        if e_s[-2] > 0 and s > 20:
            emp.append(elem)

    y = np.array(emp)
    z = np.array([y])
    return z


#cohort = "HMP1-2"
#stem = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/midas_output/"
stem = "/u/project/ngarud/wrshoema/snv_prevalence/data/midas_output_v1.2.1/"
#num = int(sys.argv[1])-1
species = sys.argv[1]
sys.stderr.write("species var is:" + str(species) + "\n")

sample = sys.argv[2]
sys.stderr.write("sample var is:" + str(sample) + "\n")

outdir = sys.argv[3]
sys.stderr.write("outdir var is:" + str(outdir) + "\n")

species = species.strip('.snps.gz')


sys.stderr.write("species var is:" + str(species) + "\n")
sys.stderr.write("sample var is:" + str(sample) + "\n")
sys.stderr.write("outdir var is:" + str(outdir) + "\n")

#####FIRST CHECK IF SPECIES/SAMPLE HAS ALREADY BEEN PROCESSED#####
filename="/u/project/ngarud/wrshoema/snv_prevalence/data/strainfinder_output/strain_dict.pickle"
with open(filename, 'rb') as handle:
    strain_dict = pickle.load(handle)
if species in strain_dict.keys():
    if sample in strain_dict[species]["North America"]:
        sys.stdout.write("passed %s for %s" % (species, sample))
        sys.exit()


#samples = os.listdir(stem)
#sample = samples[num]
#species_list = os.listdir("{stem}{sample}/snps/output/")
#sample = "SRR7658580"
#species_list = ["Prevotella_copri_61740"]

#if not os.path.isdir(f"/u/scratch/r/rwolff/strainFinder/host_sf_inputs/{cohort}/{sample}"):
#    os.makedirs(f"/u/scratch/r/rwolff/strainFinder/host_sf_inputs/{cohort}/{sample}",exist_ok=True)

# optionally, threshold on species prevalence
#species_file = pd.read_csv(f"/u/project/ngarud/Garud_lab/metagenomic_fastq_files/{cohort}/data/species/species_prevalence.txt.bz2",index_col=0,sep="\t",engine="python")
#species_file = species_file[species_file["prevalence"] >= 5]
#species_file = [e.split(".")[0] for e in os.listdir(f"{stem}{sample}/snps/output/")]
#species_list = list(species_file)
#print(species_list)

#for species in species_list:

sys.stdout.write(str(species))
snps_loc = "%s%s/snps/output/%s.snps.gz" % (stem,sample,species)
    #snps_loc = "C:/Users/sarah/Garud Lab/Prevotella_copri_61740.snps.gz"
df = pd.read_csv(snps_loc,sep='\t')
counts = np.array(df['count_atcg'])
atcg = get_counts(counts,0)

if not os.path.isdir("%s%s/" % (outdir,species)):
    os.makedirs("%s%s/" % (outdir,species))

with open("%s%s/%s.pkl" % (outdir,species,sample),'wb') as f:
    pickle.dump(atcg, f,protocol=2)

sys.stdout.write("{0} finished! \n".format(species))
sys.stdout.write("{0} finished! \n \n \n".format(sample))
sys.stdout.flush()
