import numpy
#from config import *
import parse_midas_data
import config
from collections import defaultdict
import os

###############################################################################
#
# Methods for parsing sample metadata
#
###############################################################################



flatten = lambda l: [item for sublist in l for item in sublist]


# ===========================================================================
# sample_pair_to_tp_pair: converts sample pair to general timepoint pair form
# example: frozenset('I1', 'I3')
# A for adult, M for mother, I for infant
# ===========================================================================

def sample_pair_to_tp_pair(sample_i, sample_j, sample_order_map, hmp_samples, mother_samples):
	if sample_i in hmp_samples:
		tp_i = 'A' + str(sample_order_map[sample_i][1])
	else: # If not HMP, assume mother or infant
		tp_i = ('M' if sample_i in mother_samples else 'I') + str(sample_order_map[sample_i][1])
	if sample_j in hmp_samples:
		tp_j = 'A' + str(sample_order_map[sample_j][1])
	else: # If not HMP, assume mother or infant
		tp_j = ('M' if sample_j in mother_samples else 'I') + str(sample_order_map[sample_j][1])
	tp_pair = frozenset((tp_i, tp_j))
	return tp_pair


# ===========================================================================
# load_qp_samples: returns QP status dictionary given samples, species
# 'qp', 'non-qp' [high coverage], 'low-coverage' (uses pickle)
# ===========================================================================

def load_qp_samples(desired_samples, species_name, force_repickle = False):

	import pickle, os.path
	pickle_fn = "%s/pickles/qp_samples/%s_qp_sample_dict.pkl" % (config.data_directory, species_name)

	if force_repickle or not os.path.isfile(pickle_fn):
		all_samples = get_sample_names('all')
		qp_sample_sets = calculate_qp_samples(all_samples, species_name)
		pickle.dump(qp_sample_sets, open(pickle_fn, 'wb'))
		return qp_sample_sets
	else:
		qp_sample_sets = pickle.load(open(pickle_fn, 'rb'))
		for cat in qp_sample_sets: # qp, non-qp, low-coverage
			old_sample_list = list(qp_sample_sets[cat])
			for sample in old_sample_list:
				if sample not in desired_samples:
					qp_sample_sets[cat].remove(sample)
		return qp_sample_sets

# ===========================================================================
# FUNCTIONS FOR SAMPLE-METADATA MAPS
# ===========================================================================

# ====================================================================================
# parse_sample_metadata_map
#
# Loads metadata for HMP, Yassour, Backhed, Ferretti, Shao, Olm samples
# Returns map:
# sample -> (subject_id, sample_id, accession_id, country, continent, temporal_order)
#
# By default (at least for postprocessing steps), include all samples. But can
# restrict to fecal samples only, and good timepoints only (some Shao have 'NA' tp)
# ====================================================================================

def parse_sample_metadata_map(fecal_only = False, good_tp_only = False):
	metadata_dir = config.metadata_directory
	#samples_dir = "%s/final_sample_lists" % metadata_dir
	sample_metadata_map = {} # What is returned!

	# First load HMP metadata (469 samples)

	#samples_fpath = "%s/HMP1-2_samples.txt" % samples_dir
	samples_fpath = "%s/HMP1-2_samples.txt" % metadata_dir
    #samples_fpath = "%sHMP1-2_sample_ids.txt" % metadata_dir

	hmp_samples = [line.strip() for line in open(samples_fpath, 'r')]
	hmp_samples = parse_merged_sample_names(hmp_samples) # Remove c's

	#with open("%s/HMP1-2_metadata.txt" % metadata_dir, 'r') as metadata_file:
	with open("%sHMP1-2_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			subject_id, sample_id, accession_id, country, continent, order = line.strip().split('\t')
			order = int(order)
			if sample_id in hmp_samples:
				sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, country, continent, order)

	# Then load Backhed data (391 samples)

	timept_order_map_mother = {'M':1}
	timept_order_map_infant = {'B':1,'4M':2,'12M':3}

	with open("%s/Backhed_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			accession_id, subject_id, timept = line.strip().split('\t')
			# Using the family/study_id as subject id, and specify mother (-M) vs infant (-I)
			subject_id = subject_id + ('-M' if timept == 'M' else '-I')
			sample_id = accession_id # Sample ID same as run accession
			order = timept_order_map_mother[timept] if timept == 'M' else timept_order_map_infant[timept]
			sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'Sweden', 'Europe', order)

	# Then load Ferretti data (119 fecal samples)
	timept_order_map_mother = {'t0':1}
	timept_order_map_infant = {'t1':1,'t2':2,'t3':3,'t4':4,'t5':5}

	with open("%s/Ferretti_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			items = line.strip().split("\t")
			accession_id = items[4]
			subject_id = items[6][:11] # Using first 11 characters of experiment_alias as subject id (e.g. CA_C10055IS)
			sample_id = accession_id # Sample ID same as run accession
			timept = items[5][-2:] # Using last two characters of experiment_alias to identify timepoint
			order = timept_order_map_mother[timept] if subject_id[-2:] == 'MS' else timept_order_map_infant[timept]
			if fecal_only:
				if items[6][15:17] == 'FE': # Restrict to fecal samples
					sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'Italy', 'Europe', order)
			else:
				sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'Italy', 'Europe', order)

	# Then load Yassour data (286 samples)
	timept_order_map_mother = {'Mother:Gest':1, 'Mother:Birth':2,'Mother:3 months':3}
	timept_order_map_infant = {'Child:Birth':1,'Child:14 days':2,'Child:1 month':3,'Child:2 months':4,'Child:3 months':5}

	with open("%s/Yassour_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			items = line.strip().split("\t")
			accession_id = items[4]
			subject_id = items[7][:7] # Using first 7 characters of sample_title as subject id (e.g. M0059-M)
			sample_id = accession_id # Sample ID same as run accession
			timept = items[7][6:] # Using characters after 6th of sample_title to identify timepoint
			order = timept_order_map_mother[timept] if subject_id[-1] == 'M' else timept_order_map_infant[timept]
			sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'Finland', 'Europe', order)

	# Then load Shao data (1679 samples - 1676 excluding NA timepoints)

	with open("%s/Shao_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			accession_id, _, subject_id, timept, infancy_months = line.strip().split("\t")
			sample_id = accession_id # Sample ID same as run accession

			# Adjust subject ID by appending -M or -I
			if timept == 'Mother' and infancy_months == 'Mother':
				subject_id += '-M'
				order = 1 # Only one mother timepoint
			elif timept == 'Infancy':
				subject_id += '-I'

				if infancy_months == 'NA':
					if good_tp_only == True:
						continue # Skip if infancy months field is NA
					else:
						infancy_months = -999 # Bogus negative number

				order = float(infancy_months) * 30.5 # Convert months to approx. days
			elif infancy_months == 'Neonatal':
				subject_id += '-I'
				order = int(timept) # In days since birth

			sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'United Kingdom', 'Europe', order)

	# Then load Olm data (898 samples)

	olm_samples = []
	for campaign in ['NIH1', 'NIH2', 'NIH3', 'NIH4', 'Sloan2']:
        #samples_fpath = "%s/Olm_%s_samples.txt" % (samples_dir, campaign)
		samples_fpath = "%s/Olm_%s_samples.txt" % (metadata_dir, campaign)
		olm_sub_samples = [line.strip() for line in open(samples_fpath, 'r')]
		olm_samples += list(parse_merged_sample_names(olm_sub_samples)) # Remove c's

	with open("%s/Olm_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			items = line.strip().split("\t")
			if len(items) == 10: # Must have available accession
				subject_id = items[1]
				timept = items[2]
				accession_id = items[9]
				sample_id = accession_id # Sample ID same as run accession
				order = int(timept) # In days since birth
				if accession_id in olm_samples: # Restrict to considered samples
					sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'United States', 'North America', order)

	return sample_metadata_map










# ==========================================================================
# Returns good_species_list as defined above, but from pickle file
# ==========================================================================

def load_pickled_good_species_list():

	import pickle
	pickle_path = '%s/pickles/good_species_list.txt' % config.data_directory

	if os.path.isfile(pickle_path):
		return pickle.load(open(pickle_path, 'rb'))
	else:
		good_species_list = parse_good_species_list()
		pickle.dump(good_species_list, open(pickle_path, 'wb'))
		return good_species_list



# ===========================================================================
# FUNCTIONS FOR DEALING WITH SAMPLE LISTS
# ===========================================================================

# ===========================================================================
# get_sample_names: returns list of samples in specified cohort/timepoint
#
# cohort: one of 'HMP', 'Backhed', 'Ferretti', 'Yassour', Shao', 'Olm'
# (case insensitive),	'mother', 'infant', 'all', or 'all-dict'
# note that 'dict' is (cohort->tp->samples)
#
# timepoint:
# 	HMP - 1, 2, 3
#		Backhed - B, 4M, 12M, M
#		Ferretti - M0, I1, I2, I3, I4, I5
#		Yassour - MGest, MBirth, M3, CBirth, C14, C1, C2, C3
# 	Shao - Mother, Neonatal (4-21 days), Infancy
# 	Olm - NIH1, NIH2, NIH3, NIH4, Sloan2
#		or 'all' for every sample in a cohort
# ===========================================================================

def get_sample_names(cohort, timepoint = 'all'):

	sample_dict = {'hmp': defaultdict(set), 'backhed': defaultdict(set), \
								 'ferretti': defaultdict(set), 'yassour': defaultdict(set), \
								 'shao': defaultdict(set), 'olm': defaultdict(set)}

	# Note: we need the sample lists because not all samples in the metadata
	# files should be returned. The final ID is based on metadata file, however
	# (so -c suffixes for combined samples will be omitted.)

	metadata_dir = config.metadata_directory
	samples_dir = "%sfinal_sample_lists" % metadata_dir

	# HMP - expect 469 samples
	#samples_fpath = "%sHMP1-2_samples.txt" % samples_dir
	samples_fpath = "%sHMP1-2_samples.txt" % metadata_dir
    #samples_fpath = "%sHMP1-2_sample_ids.txt" % metadata_dir

	hmp_samples = [line.strip() for line in open(samples_fpath, 'r')]
	hmp_samples = parse_merged_sample_names(hmp_samples) # Remove c's

	metadata_fpath = "%s/HMP1-2_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]

	for sample in metadata[1:]:
		_, sample_id, _, _, _, tp = sample
		if sample_id in hmp_samples:
			sample_dict['hmp'][tp].add(sample_id)

	# Backhed - expect 391 samples
	#samples_fpath = "%sBackhed_samples.txt" % samples_dir
	samples_fpath = "%sBackhed_samples.txt" % metadata_dir
	backhed_samples = [line.strip() for line in open(samples_fpath, 'r')]
	backhed_samples = parse_merged_sample_names(backhed_samples) # Remove c's

	metadata_fpath = "%sBackhed_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]

	for sample in metadata[1:]:
		sample_id, _, tp = sample
		if sample_id in backhed_samples:
			sample_dict['backhed'][tp].add(sample_id)

	# Ferretti - expect 119 samples
	samples_fpath = "%sFerretti_samples.txt" % metadata_dir
	ferretti_samples = [line.strip() for line in open(samples_fpath, 'r')]
	ferretti_samples = parse_merged_sample_names(ferretti_samples) # Remove c's

	metadata_fpath = "%sFerretti_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]

	for sample in metadata[1:]:
		sample_id = sample[4]
		tp = sample[6][9] + sample[6][19]
		# Restrict to fecal samples
		body_site = sample[6][15:17]
		if body_site == 'FE' and sample_id in ferretti_samples:
			sample_dict['ferretti'][tp].add(sample_id)

	# Yassour - expect 286 samples
	samples_fpath = "%sYassour_samples.txt" % metadata_dir
	yassour_samples = [line.strip() for line in open(samples_fpath, 'r')]
	yassour_samples = parse_merged_sample_names(yassour_samples) # Remove c's

	metadata_fpath = "%sYassour_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]

	for sample in metadata[1:]:
		sample_id = sample[4]
		tp_raw = sample[7][6:].split(':')
		tp = tp_raw[0][0] + tp_raw[1].split()[0]
		if sample_id in yassour_samples:
			sample_dict['yassour'][tp].add(sample_id)

	# Shao - expect 1679 samples
	samples_fpath = "%sShao_samples.txt" % metadata_dir
	shao_samples = [line.strip() for line in open(samples_fpath, 'r')]
	shao_samples = parse_merged_sample_names(shao_samples) # Remove c's

	metadata_fpath = "%sShao_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]

	for sample in metadata[1:]:
		sample_id, _, _, neonatal_tp, tp_cat = sample
		if tp_cat == 'NA':
			continue # Ignore if order info not known
		tp_cat = 'Infancy' if neonatal_tp == 'Infancy' else tp_cat
		if sample_id in shao_samples:
			sample_dict['shao'][tp_cat].add(sample_id)

	# Olm - expect 898 samples
	olm_samples = []
	for campaign in ['NIH1', 'NIH2', 'NIH3', 'NIH4', 'Sloan2']:
        #samples_fpath = "%sOlm_%s_samples.txt" % (samples_dir, campaign)
		samples_fpath = "%sOlm_%s_samples.txt" % (metadata_dir, campaign)
		olm_sub_samples = [line.strip() for line in open(samples_fpath, 'r')]
		olm_samples += list(parse_merged_sample_names(olm_sub_samples)) # Remove c's

	metadata_fpath = "%sOlm_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]

	for sample in metadata[1:]:
		try:
			sample_id = sample[9]
			campaign = sample[3]
			if sample_id in olm_samples:
				sample_dict['olm'][campaign].add(sample_id)
		except:
			continue

	# Aggregate + convert sets to lists
	all_samples = []
	for lcohort in sample_dict:
		for ltp in sample_dict[lcohort]:
			sample_dict[lcohort][ltp] = list(sample_dict[lcohort][ltp])
			all_samples += sample_dict[lcohort][ltp]

	# All HMP samples
	hmp_samples = flatten([sample_dict['hmp'][htp] for htp in sample_dict['hmp']])

	# All mother and infant samples
	mother_samples = []
	mother_samples += sample_dict['backhed']['M']
	mother_samples += sample_dict['ferretti']['M0']
	mother_samples += sample_dict['shao']['Mother']
	mother_samples += flatten([sample_dict['yassour'][ytp] for ytp in ['MGest', 'MBirth', 'M3']])
	infant_samples = [s for s in all_samples if (s not in mother_samples and s not in hmp_samples)]

	general_cohort_dict = {'all': all_samples, 'mother': mother_samples, 'infant': infant_samples}

	cohort = cohort.lower() # for case insensitivity

	if cohort in general_cohort_dict: # all, mother, infant
		return general_cohort_dict[cohort]
	elif cohort == 'all-dict': # all-dict
		return sample_dict
	elif timepoint == 'all': # specific cohort, all
		all_cohort_samples = []
		for ltp in sample_dict[cohort]:
			all_cohort_samples += sample_dict[cohort][ltp]
		return all_cohort_samples
	else: # specific cohort, specific timepoint
		return sample_dict[cohort][timepoint]





# ====================================================================================
# extract_sample_metadata_map
#
# Using passed in sample-metadata map, extracts information for only one column
# (options: subject_id, country, continent or order) and returns new map
# Loads the default full sample-metadata map if nothing passed in
# ====================================================================================

def extract_sample_metadata_map(field, sample_metadata_map = None):

	field_dict = {"subject_id": 0, "country": 3, "continent": 4, "order": 5}
	field_idx = field_dict[field] if (field in field_dict) else 0 # Defaults to subject_id

	if sample_metadata_map is None:
		sample_metadata_map = parse_sample_metadata_map() # Load it

	extracted_sample_metadata_map = {}
	for sample in sample_metadata_map:
		extracted_sample_metadata_map[sample] = sample_metadata_map[sample][field_idx]

	return extracted_sample_metadata_map



# ====================================================================================
# parse_sample_subject_map, parse_sample_country_map, parse_sample_continent_map
#
# Convenience functions for extract_sample_metadata_map
# ====================================================================================

def parse_sample_subject_map(sample_metadata_map = None):
	return extract_sample_metadata_map("subject_id", sample_metadata_map)

def parse_sample_country_map(sample_metadata_map = None):
	return extract_sample_metadata_map("country", sample_metadata_map)

def parse_sample_continent_map(sample_metadata_map = None):
	return extract_sample_metadata_map("continent", sample_metadata_map)








# ====================================================================================
# parse_sample_order_map
#
# Returns map from sample -> (subject_id, temporal_order)
# ====================================================================================

def parse_sample_order_map(sample_metadata_map = None):

	if sample_metadata_map is None:
		sample_metadata_map = parse_sample_metadata_map() # Load it

	sample_order_map = {}
	for sample in sample_metadata_map:
			subject_id, _, _, _, _, order = sample_metadata_map[sample]
			sample_order_map[sample] = (subject_id, order)

	return sample_order_map


def parse_merged_sample_names(items):
    samples = []
    for item in items:
        sample = item.strip()
        if sample.endswith('c'):
            sample = sample[:-1]
        samples.append(sample)

    samples = numpy.array(samples)
    return samples


def calculate_sample_subject_map(subject_sample_map):

    # invert subject sample map
    sample_subject_map = {}
    for subject in subject_sample_map.keys():
        for sample in subject_sample_map[subject].keys():
            sample_subject_map[sample] = subject

    return sample_subject_map


###############################################################################
#
# Creates a map of indexes from one list of samples (sample_list_from)
# to another list of samples (sample_list_to). The from list must be a
# strict subset of the to list.
#
###############################################################################
def calculate_sample_idx_map(sample_list_from, sample_list_to):

    sample_list_to = list(sample_list_to)
    sample_map = {}
    for i in xrange(0,len(sample_list_from)):
        sample_map[i] = sample_list_to.index(sample_list_from[i])

    return sample_map

def apply_sample_index_map_to_indices(sample_idx_map, idxs):
    new_idxs = (numpy.array([sample_idx_map[i] for i in idxs[0]]), numpy.array([sample_idx_map[i] for i in idxs[1]]))
    return new_idxs

def sample_name_lookup(sample_name, samples):

    for sample in samples:

        if sample.startswith(sample_name):
            return sample

    return ""

###############################################################################
#
# Prunes sample list to remove multiple timepoints from same subject
# Returns len(sampe_list) boolean array with element=False if sample was pruned
#
###############################################################################
def calculate_unique_samples(subject_sample_map, sample_list=[]):

    if len(sample_list)==0:
        sample_list = list(sorted(flatten_samples(subject_sample_map).keys()))

    # invert subject sample map
    sample_subject_map = {}
    for subject in subject_sample_map.keys():
        for sample in subject_sample_map[subject].keys():
            sample_subject_map[sample] = subject

    subject_idx_map = {}

    for i in xrange(0,len(sample_list)):
        sample = sample_list[i]
        if sample.endswith('c'):
            sample = sample[:-1]
        subject = sample_subject_map[sample]
        if not subject in subject_idx_map:
            subject_idx_map[subject] = i

    unique_idxs = numpy.zeros(len(sample_list),dtype=numpy.bool_)
    for i in subject_idx_map.values():
        unique_idxs[i]=True

    return unique_idxs

###
#
# Returns a vector of true false values
#
###
def calculate_samples_in_different_subjects(subject_sample_map, sample_list, focal_sample):

    # invert subject sample map
    sample_subject_map = {}
    for subject in subject_sample_map.keys():
        for sample in subject_sample_map[subject].keys():
            sample_subject_map[sample] = subject

    in_different_subject = []
    for sample in sample_list:
        if sample_subject_map[sample] == sample_subject_map[focal_sample]:
            in_different_subject.append(False)
        else:
            in_different_subject.append(True)

    return numpy.array(in_different_subject)


###############################################################################
#
# Prunes sample list to only include samples from allowed countries
# Returns len(sampe_list) boolean array with element=False if sample was pruned
#
###############################################################################
def calculate_country_samples(sample_country_map, sample_list=[], allowed_countries=set([])):

    if len(sample_list)==0:
        sample_list = list(sorted(sample_country_map.keys()))

    allowed_idxs = []
    for sample in sample_list:

        if sample.endswith('c'):
            desired_sample = sample[:-1]
        else:
            desired_sample = sample

        if (len(allowed_countries))==0 or (sample_country_map[desired_sample] in allowed_countries):
            allowed_idxs.append(True)
        else:
            allowed_idxs.append(False)

    allowed_idxs = numpy.array(allowed_idxs)
    return allowed_idxs


###############################################################################
#
# For a given list of samples, calculates which belong to different subjects
# which belong to different timepoints in same subject, and which are the same
# timepoint.
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs,
# each of which is a tuple with idx1 and idx2. All pairs are included
# only once.
#
###############################################################################
def calculate_subject_pairs(subject_sample_map, sample_list=[]):

    if len(sample_list)==0:
        sample_list = list(sorted(flatten_samples(subject_sample_map).keys()))

    new_sample_list = []
    for sample in sample_list:
        if sample.endswith('c'):
            new_sample_list.append(sample[:-1])
        else:
            new_sample_list.append(sample)

    sample_list = new_sample_list

    # invert subject sample map
    sample_subject_map = {}
    for subject in subject_sample_map.keys():
        for sample in subject_sample_map[subject].keys():
            sample_subject_map[sample] = subject

    same_sample_idx_lower = []
    same_sample_idx_upper = []
    same_subject_idx_lower = []
    same_subject_idx_upper = []
    diff_subject_idx_lower = []
    diff_subject_idx_upper = []

    for i in xrange(0,len(sample_list)):
        same_sample_idx_lower.append(i)
        same_sample_idx_upper.append(i)
        for j in xrange(0,i):
            if sample_subject_map[sample_list[i]]==sample_subject_map[sample_list[j]]:
                same_subject_idx_lower.append(i)
                same_subject_idx_upper.append(j)
            else:
                diff_subject_idx_lower.append(i)
                diff_subject_idx_upper.append(j)

    same_sample_idxs = (numpy.array(same_sample_idx_lower,dtype=numpy.int32), numpy.array(same_sample_idx_upper,dtype=numpy.int32))

    same_subject_idxs = (numpy.array(same_subject_idx_lower,dtype=numpy.int32), numpy.array(same_subject_idx_upper,dtype=numpy.int32))

    diff_subject_idxs = (numpy.array(diff_subject_idx_lower,dtype=numpy.int32), numpy.array(diff_subject_idx_upper,dtype=numpy.int32))

    return same_sample_idxs, same_subject_idxs, diff_subject_idxs

###############################################################################
#
# For a given list of samples, calculates which belong to different subjects
# which belong to different timepoints in same subject, and which are the same
# timepoint.
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs,
# each of which is a tuple with idx1 and idx2. All pairs are included
# only once.
#
###############################################################################
def calculate_old_ordered_subject_pairs(sample_order_map, sample_list=[]):

    same_sample_idx_lower = []
    same_sample_idx_upper = []
    same_subject_idx_lower = []
    same_subject_idx_upper = []
    diff_subject_idx_lower = []
    diff_subject_idx_upper = []

    for i in xrange(0,len(sample_list)):
        for j in xrange(i,len(sample_list)):
            # loop over all pairs of samples

            if i==j:
                same_sample_idx_lower.append(i)
                same_sample_idx_upper.append(j)
            else:

                subject1, order1 = sample_order_map[sample_list[i]]
                subject2, order2 = sample_order_map[sample_list[j]]

                if subject1==subject2:
                    # same subject!
                    if order2-order1==1:
                        # consecutive samples
                        same_subject_idx_lower.append(i)
                        same_subject_idx_upper.append(j)
                    elif order1-order2==1:
                        # consecutive samples
                        same_subject_idx_lower.append(j)
                        same_subject_idx_upper.append(i)
                    else:
                        # do not add
                        pass

                else:
                    # different subjects!
                    # Only take first one (to prevent multiple comparisons)
                    if order1==1 and order2==1:
                        diff_subject_idx_lower.append(i)
                        diff_subject_idx_upper.append(j)

    same_sample_idxs = (numpy.array(same_sample_idx_lower,dtype=numpy.int32), numpy.array(same_sample_idx_upper,dtype=numpy.int32))

    same_subject_idxs = (numpy.array(same_subject_idx_lower,dtype=numpy.int32), numpy.array(same_subject_idx_upper,dtype=numpy.int32))

    diff_subject_idxs = (numpy.array(diff_subject_idx_lower,dtype=numpy.int32), numpy.array(diff_subject_idx_upper,dtype=numpy.int32))

    return same_sample_idxs, same_subject_idxs, diff_subject_idxs

###############################################################################
#
# For a given list of samples, calculates which belong to different subjects
# which belong to different timepoints in same subject, and which are the same
# timepoint.
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs,
# each of which is a tuple with idx1 and idx2. All pairs are included
# only once.
#
###############################################################################
def calculate_ordered_subject_pairs(sample_order_map, sample_list=[], within_host_type='consecutive'):

    same_sample_idx_lower = []
    same_sample_idx_upper = []
    same_subject_idx_lower = []
    same_subject_idx_upper = []
    diff_subject_idx_lower = []
    diff_subject_idx_upper = []

    diff_subject_pair_map = {}
    same_subject_pair_map = {}
    # in this list, store lower_idx, upper_idx, interval_length

    # same sample pairs.. trivial
    same_sample_idx_lower = numpy.arange(0,len(sample_list))
    same_sample_idx_upper = numpy.arange(0,len(sample_list))

    # reconstruct "timeseries" for each subject
    subject_order_idx_map = {}
    for i in xrange(0,len(sample_list)):
        subject, order = sample_order_map[sample_list[i]]

        if subject not in subject_order_idx_map:
            subject_order_idx_map[subject] = {}

        subject_order_idx_map[subject][order] = i

    # create index pairs within subjects
    for subject in subject_order_idx_map:

        sorted_orders = list(sorted(subject_order_idx_map[subject].keys()))

        # if at least two samples
        if len(sorted_orders)>1.5:

            if within_host_type=='longest':
                same_subject_idx_lower.append( subject_order_idx_map[subject][sorted_orders[0]] )
                same_subject_idx_upper.append( subject_order_idx_map[subject][sorted_orders[-1]] )
            elif within_host_type=='consecutive':
                for order_idx in xrange(1,len(sorted_orders)):
                    same_subject_idx_lower.append( subject_order_idx_map[subject][sorted_orders[order_idx-1]] )
                    same_subject_idx_upper.append( subject_order_idx_map[subject][sorted_orders[order_idx]] )
            elif within_host_type=='nonconsecutive':
                for order_idx_i in xrange(0,len(sorted_orders)):
                    for order_idx_j in xrange(order_idx_i+1,len(sorted_orders)):
                        same_subject_idx_lower.append( subject_order_idx_map[subject][sorted_orders[order_idx_i]] )
                        same_subject_idx_upper.append( subject_order_idx_map[subject][sorted_orders[order_idx_j]] )

    # now create index pairs in different subjects
    sorted_subjects = sorted(subject_order_idx_map.keys())

    for subject_i_idx in xrange(0,len(sorted_subjects)):
        subject_i = sorted_subjects[subject_i_idx]
        i = subject_order_idx_map[subject_i][min(subject_order_idx_map[subject_i].keys())]

        for subject_j_idx in xrange(subject_i_idx+1,len(sorted_subjects)):
            subject_j = sorted_subjects[subject_j_idx]
            j = subject_order_idx_map[subject_j][min(subject_order_idx_map[subject_j].keys())]

            diff_subject_idx_lower.append(i)
            diff_subject_idx_upper.append(j)

    same_sample_idxs = (numpy.array(same_sample_idx_lower,dtype=numpy.int32), numpy.array(same_sample_idx_upper,dtype=numpy.int32))

    same_subject_idxs = (numpy.array(same_subject_idx_lower,dtype=numpy.int32), numpy.array(same_subject_idx_upper,dtype=numpy.int32))

    diff_subject_idxs = (numpy.array(diff_subject_idx_lower,dtype=numpy.int32), numpy.array(diff_subject_idx_upper,dtype=numpy.int32))

    return same_sample_idxs, same_subject_idxs, diff_subject_idxs


###############################################################################
#
# For a given list of samples, calculates which belong to different subjects
# which belong to different timepoints in same subject, and which are the same
# timepoint.
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs,
# each of which is a tuple with idx1 and idx2.
# The same pair can be included multiple times if there are multiple timepoints
#
###############################################################################
def calculate_nonconsecutive_ordered_subject_pairs(sample_order_map, sample_list=[]):

    same_sample_idx_lower = []
    same_sample_idx_upper = []
    same_subject_idx_lower = []
    same_subject_idx_upper = []
    diff_subject_idx_lower = []
    diff_subject_idx_upper = []

    for i in xrange(0,len(sample_list)):
        for j in xrange(i,len(sample_list)):
            # loop over all pairs of samples

            if i==j:
                same_sample_idx_lower.append(i)
                same_sample_idx_upper.append(j)
            else:

                subject1, order1 = sample_order_map[sample_list[i]]
                subject2, order2 = sample_order_map[sample_list[j]]

                if subject1==subject2:
                    # same subject!
                    if order2-order1>0.5:
                        # consecutive samples
                        same_subject_idx_lower.append(i)
                        same_subject_idx_upper.append(j)
                    elif order1-order2>0.5:
                        # consecutive samples
                        same_subject_idx_lower.append(j)
                        same_subject_idx_upper.append(i)
                    else:
                        # do not add
                        pass

                else:
                    # different subjects!
                    # Only take first one (to prevent multiple comparisons)
                    if order1==1 and order2==1:
                        diff_subject_idx_lower.append(i)
                        diff_subject_idx_upper.append(j)

    same_sample_idxs = (numpy.array(same_sample_idx_lower,dtype=numpy.int32), numpy.array(same_sample_idx_upper,dtype=numpy.int32))

    same_subject_idxs = (numpy.array(same_subject_idx_lower,dtype=numpy.int32), numpy.array(same_subject_idx_upper,dtype=numpy.int32))

    diff_subject_idxs = (numpy.array(diff_subject_idx_lower,dtype=numpy.int32), numpy.array(diff_subject_idx_upper,dtype=numpy.int32))

    return same_sample_idxs, same_subject_idxs, diff_subject_idxs


###############################################################################
# For a given list of samples, calculates which belong to different timepoints in same subject
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs,
# each of which is a tuple with idx1 and idx2. All pairs are included
# only once.
#
###############################################################################
def calculate_ordered_subject_triplets(sample_order_map, sample_list=[]):

    same_subject_idxs = []

    for i in xrange(0,len(sample_list)):

        subject1, order1 = sample_order_map[sample_list[i]]
        if order1 != 1:
            continue

        for j in xrange(0,len(sample_list)):

            subject2, order2 = sample_order_map[sample_list[j]]

            if subject2 != subject1:
                continue

            if order2 != 2:
                continue

            for k in xrange(0,len(sample_list)):

                subject3, order3 = sample_order_map[sample_list[k]]

                if subject3 != subject1:
                    continue

                if order3 != 3:
                    continue

                # if you get here, a triplet!
                same_subject_idxs.append((i,j,k))

    return same_subject_idxs



###############################################################################
#
# Calculates the subset that are sampled three times
#
###############################################################################
def calculate_triple_samples(sample_order_map, sample_list=[]):

    sample_idx_map = {}


    for i in xrange(0,len(sample_list)):
        subject, order = sample_order_map[sample_list[i]]

        if subject not in sample_idx_map:
            sample_idx_map[subject] = {}

        sample_idx_map[subject][order] = i

    triple_samples = []
    for subject in sample_idx_map.keys():
        if len(sample_idx_map[subject].keys()) > 2:
            triple_samples.append(numpy.array([sample_idx_map[subject][order] for order in sorted(sample_idx_map[subject].keys())]))

    return triple_samples


###############################################################################
#
# Returns matrix where rows are subjects and columns are hosts
# where A_ih = 1 if sample i is in host h
#
###############################################################################
def calculate_sample_subject_matrix(samples):

    sample_idx_map = {samples[i]:i for i in xrange(0,len(samples))}

    subject_sample_map = parse_subject_sample_map()
    subjects = subject_sample_map.keys()

    sample_subject_matrix = numpy.zeros((len(samples),len(subjects)),dtype=numpy.bool)

    for subject_idx in xrange(0,len(subjects)):
        for sample in subject_sample_map[subjects[subject_idx]]:
            if sample in sample_idx_map:
                sample_subject_matrix[sample_idx_map[sample], subject_idx] = True

    return sample_subject_matrix, subjects





    same_sample_idx_lower = []
    same_sample_idx_upper = []
    same_subject_idx_lower = []
    same_subject_idx_upper = []
    diff_subject_idx_lower = []
    diff_subject_idx_upper = []

    for i in xrange(0,len(sample_list)):
        for j in xrange(i,len(sample_list)):
            # loop over all pairs of samples

            if i==j:
                same_sample_idx_lower.append(i)
                same_sample_idx_upper.append(j)
            else:

                subject1, order1 = sample_order_map[sample_list[i]]
                subject2, order2 = sample_order_map[sample_list[j]]

                if subject1==subject2:
                    # same subject!
                    if order2-order1==1:
                        # consecutive samples
                        same_subject_idx_lower.append(i)
                        same_subject_idx_upper.append(j)
                    elif order1-order2==1:
                        # consecutive samples
                        same_subject_idx_lower.append(j)
                        same_subject_idx_upper.append(i)
                    else:
                        # do not add
                        pass

                else:
                    # different subjects!
                    # Only take first one (to prevent multiple comparisons)
                    if order1==1 and order2==1:
                        diff_subject_idx_lower.append(i)
                        diff_subject_idx_upper.append(j)

    same_sample_idxs = (numpy.array(same_sample_idx_lower,dtype=numpy.int32), numpy.array(same_sample_idx_upper,dtype=numpy.int32))

    same_subject_idxs = (numpy.array(same_subject_idx_lower,dtype=numpy.int32), numpy.array(same_subject_idx_upper,dtype=numpy.int32))

    diff_subject_idxs = (numpy.array(diff_subject_idx_lower,dtype=numpy.int32), numpy.array(diff_subject_idx_upper,dtype=numpy.int32))

    return same_sample_idxs, same_subject_idxs, diff_subject_idxs

###############################################################################
#
# Returns a flat map of all the replicate sets for
# the samples in subject_sample_map, indexed by sample key
#
###############################################################################
def flatten_samples(subject_sample_map):

    grouping_replicate_map = {}
    for subject in sorted(subject_sample_map.keys()):
        for sample in sorted(subject_sample_map[subject].keys()):
            grouping_replicate_map[sample] = subject_sample_map[subject][sample]

    return grouping_replicate_map


###############################################################################
#
# Returns a flat map of the merged replicate sets for each subject,
# indexed by subject key
#
###############################################################################
def flatten_subjects(subject_sample_map):

    grouping_replicate_map = {}
    for subject in sorted(subject_sample_map.keys()):
        merged_replicates = set()
        for sample in subject_sample_map[subject].keys():
            merged_replicates.update(subject_sample_map[subject][sample])
        grouping_replicate_map[subject] = merged_replicates

    return grouping_replicate_map


###############################################################################
#
# groupings = ordered list of nonoverlapping sets of sample names
# samples = ordered list of samples
#
# returns: list whose i-th element contains a numpy array of idxs
#          of the items in samples that are present in the ith grouping
#
###############################################################################
def calculate_grouping_idxs(groupings, samples):

    grouping_idxs = []
    for i in xrange(0,len(groupings)):

        idxs = []
        for j in xrange(0,len(samples)):
            if samples[j] in groupings[i]:
                idxs.append(j)
        idxs = numpy.array(idxs,dtype=numpy.int32)
        #print idxs
        grouping_idxs.append(idxs)

    return grouping_idxs
