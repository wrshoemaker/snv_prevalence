from __future__ import division
import os
import sys


#species = sys.argv[1]

species = 'Bacteroides_vulgatus_57955'


dir = "/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/midas_output_bam/%s/" % species


#dir = "/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/%s/" % 'Bacteroides_vulgatus_57955'


# loop through BED file
bed_filepath = "/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/rep_genomes/%s/genome.bed" % species
file_bed = open(bed_filepath, "r")
scaffolds_to_keep = []
#file_bed.readline() # header
for line in file_bed:
    line_split = line.strip().split('\t')
    scaffolds_to_keep.append(line_split[0])

file_bed.close()

print(scaffolds_to_keep)

# delete all empty files

for filename in os.listdir(dir):

    filepath = '%s%s' % (dir, filename)
    size = os.path.getsize(filepath)

    if size == 0:

        os.remove(filepath)



# now

for filename in os.listdir(dir):

    if filename.endswith(".header"):
        # print(os.path.join(directory, filename))
        file = open("%s%s" % (dir, filename), "r")
        first_line = file.readline() # header
        #for line in file:

        file.close()

    else:
        continue


print(first_line)
