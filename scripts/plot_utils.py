import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import parse_midas_data as pmd
import midas_db_utils
from collections import defaultdict

import matplotlib.cm as cm
from matplotlib import colors

good_species = 'Eubacterium_rectale_56927'
bad_species = 'Bacteroides_vulgatus_57955'

good_bad_color_dict = {good_species: 'dodgerblue', bad_species: '#FF6347'}
good_bad_color_map_dict = {good_species: cm.Blues, bad_species: cm.Reds}

variant_color_dict = {'4D': 'dodgerblue', '1D': '#FF6347'}
variant_cmap_dict = {'4D': 'Blues', '1D': 'Reds'}
sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i']

n_points=1000
color_radius=2


species_color_map_genus = {'Alistipes_finegoldii_56071': colors.hex2color(colors.cnames['lightgreen']),
                    'Alistipes_onderdonkii_55464': colors.hex2color(colors.cnames['mediumseagreen']),
                    'Alistipes_putredinis_61533': colors.hex2color(colors.cnames['seagreen']),
                    'Alistipes_shahii_62199':  colors.hex2color(colors.cnames['darkgreen']),
                    'Bacteroidales_bacterium_58650': colors.hex2color(colors.cnames['khaki']),
                    'Bacteroides_caccae_53434': colors.hex2color(colors.cnames['mediumblue']),
                    'Bacteroides_cellulosilyticus_58046': colors.hex2color(colors.cnames['aquamarine']),
                    'Bacteroides_fragilis_54507': colors.hex2color(colors.cnames['steelblue']),
                    'Bacteroides_ovatus_58035': colors.hex2color(colors.cnames['cyan']),
                    'Bacteroides_stercoris_56735': colors.hex2color(colors.cnames['lightskyblue']),
                    'Bacteroides_thetaiotaomicron_56941': colors.hex2color(colors.cnames['darkturquoise']),
                    'Bacteroides_uniformis_57318': colors.hex2color(colors.cnames['teal']),
                    'Bacteroides_vulgatus_57955': colors.hex2color(colors.cnames['dodgerblue']),
                    'Bacteroides_xylanisolvens_57185': colors.hex2color(colors.cnames['cadetblue']),
                    'Barnesiella_intestinihominis_62208': colors.hex2color(colors.cnames['darkorange']),
                    'Dialister_invisus_61905': colors.hex2color(colors.cnames['sandybrown']),
                    'Eubacterium_rectale_56927': colors.hex2color(colors.cnames['olive']),
                    'Oscillibacter_sp_60799': colors.hex2color(colors.cnames['saddlebrown']),
                    'Parabacteroides_distasonis_56985': colors.hex2color(colors.cnames['orchid']),
                    'Parabacteroides_merdae_56972': colors.hex2color(colors.cnames['darkmagenta']),
                    'Ruminococcus_bicirculans_59300': colors.hex2color(colors.cnames['orangered']),
                    'Ruminococcus_bromii_62047': colors.hex2color(colors.cnames['darkred'])}





latex_species_name_dict = {'Faecalibacterium_prausnitzii_61481':r'$Faecalibacterium \; prausnitzii$',
							'Faecalibacterium_cf_62236': r'$Faecalibacterium \, cf.$',
							'Bacteroides_salyersiae_54873': r'$Bacteroides\, salyersiae$',
							'Bacteroides_intestinalis_61596': r'$Bacteroides \, intestinalis$',
							'Bacteroides_uniformis_57318': r'$Bacteroides \, uniformis$',
							'Bacteroides_coprocola_61586': r'$Bacteroides \, coprocola$',
							'Paraprevotella_clara_33712': r'$Paraprevotella \, clara$',
							'Phascolarctobacterium_sp_59817': r'$Phascolarctobacterium \, sp.$',
							'Sutterella_wadsworthensis_56828': r'$Sutterella \, wadsworthensis$',
							'Bacteroides_fragilis_54507': r'$Bacteroides \, fragilis$',
							'Prevotella_copri_61740': r'$Prevotella \, copri$',
							'Butyrivibrio_crossotus_61674': r'$Butyrivibrio \, crossotus$',
							'Bacteroides_thetaiotaomicron_56941': r'$Bacteroides \, thetaiotaomicron$',
							'Bacteroides_ovatus_58035': r'$Bacteroides \, ovatus$',
							'Bacteroides_massiliensis_44749': r'$Bacteroides \, massiliensis$',
							'Faecalibacterium_prausnitzii_57453': r'$Faecalibacterium \, prausnitzii$',
							'Bacteroides_cellulosilyticus_58046': r'$Bacteroides \, cellulosilyticus$',
							'Oscillibacter_sp_60799': r'$Oscillibacter \, sp.$',
							'Bacteroides_finegoldii_57739': r'$Bacteroides \, finegoldii$',
							'Clostridium_sp_61482': r'$Clostridium \, sp.$',
							'Parabacteroides_merdae_56972': r'$Parabacteroides \, merdae$',
							'Odoribacter_splanchnicus_62174': r'$Odoribacter \, splanchnicus$',
							'Eubacterium_eligens_61678': r'$Eubacterium \, eligens$',
							'Alistipes_finegoldii_56071': r'$Alistipes \, finegoldii$',
							'Akkermansia_muciniphila_55290': r'$Akkermansia \, muciniphila$',
							'Roseburia_inulinivorans_61943': r'$Roseburia \, inulinivorans$',
							'Alistipes_sp_60764': r'$Alistipes \, sp.$',
							'Roseburia_intestinalis_56239': r'$Roseburia \, intestinalis$',
							'Bacteroides_xylanisolvens_57185': r'$Bacteroides \, xylanisolvens$',
							'Lachnospiraceae_bacterium_51870': r'$Lachnospiraceae \, bacterium$',
							'Parabacteroides_johnsonii_55217': r'$Parabacteroides \, johnsonii$',
							'Eubacterium_siraeum_57634': r'$Eubacterium \, siraeum$',
							'Alistipes_onderdonkii_55464': r'$Alistipes \, onderdonkii$',
							'Bacteroides_stercoris_56735': r'$Bacteroides \, stercoris$',
							'Dialister_invisus_61905': r'$Dialister \, invisus$',
							'Eubacterium_rectale_56927': r'$Eubacterium \, rectale$',
							'Sutterella_wadsworthensis_62218': r'$Sutterella \, wadsworthensis$',
							'Bacteroidales_bacterium_58650': r'$Bacteroidales \, bacterium$',
							'Barnesiella_intestinihominis_62208': r'$Barnesiella \, intestinihominis$',
							'Bacteroides_eggerthii_54457': r'$Bacteroides \, eggerthii$',
							'Bacteroides_caccae_53434': r'$Bacteroides \, caccae$',
							'Alistipes_putredinis_61533': r'$Alistipes \, putredinis$',
							'Alistipes_shahii_62199': r'$Alistipes \, shahii$',
							'Bacteroides_plebeius_61623': r'$Bacteroides \, plebeius$',
							'Butyricimonas_virosa_58742': r'$Butyricimonas \, virosa$',
							'Parabacteroides_distasonis_56985': r'$Parabacteroides \, distasonis$',
							'Ruminococcus_bicirculans_59300': r'$Ruminococcus \, bicirculans$',
							'Bacteroides_faecis_58503': r'$Bacteroides \, faecis$',
							'Bacteroides_vulgatus_57955': r'$Bacteroides \, vulgatus$'}





def get_species_color_map(all_species = pmd.parse_species_list()):

	gfo_phylum_map = midas_db_utils.load_gfo_phylum_map()

	all_phyla = set()
	phylum_species_map = defaultdict(list)

	for species in all_species:
		gfo = species.split('_')[0] # Either genus, family or order
		try:
			phylum = gfo_phylum_map[gfo]
			all_phyla.add(phylum)
			phylum_species_map[phylum].append(species)
		except:
			continue
			#print(species) # Guyana massiliensis is unclassified, skip for now

	# There are the following 7 phyla in these gut microbiota
	species_color_map = {}
	ordered_species_list = []

	# Firmicutes (81)
	species_list = sorted(phylum_species_map['Firmicutes'], key=lambda x: x[1])
	ordered_species_list += species_list
	m = get_cm_ScalerMappable(matplotlib.cm.PuRd, len(species_list), 7, 3)
	for i in range(len(species_list)):
		species_color_map[species_list[i]] = m.to_rgba(i)

	# Bacteroidetes (54)
	species_list = sorted(phylum_species_map['Bacteroidetes'], key=lambda x: x[1])
	ordered_species_list += species_list
	m = get_cm_ScalerMappable(matplotlib.cm.YlGnBu, len(species_list), 5, 2)
	for i in range(len(species_list)):
		species_color_map[species_list[i]] = m.to_rgba(i)

	# Proteobacteria (10)
	species_list = sorted(phylum_species_map['Proteobacteria'], key=lambda x: x[1])
	ordered_species_list += species_list
	m = get_cm_ScalerMappable(matplotlib.cm.Oranges, len(species_list), 2, 3)
	for i in range(len(species_list)):
		species_color_map[species_list[i]] = m.to_rgba(i)

	# Actinobacteria (9)
	species_list = sorted(phylum_species_map['Actinobacteria'], key=lambda x: x[1])
	ordered_species_list += species_list
	m = get_cm_ScalerMappable(matplotlib.cm.Wistia, len(species_list), 2, 4)
	for i in range(len(species_list)):
		species_color_map[species_list[i]] = m.to_rgba(i)

	# Fusobacteria (1), Spirochaetes (1), Verrucomicrobia (1)
	species_list = phylum_species_map['Fusobacteria'] + phylum_species_map['Spirochaetes'] + phylum_species_map['Verrucomicrobia']
	ordered_species_list += species_list
	m = get_cm_ScalerMappable(matplotlib.cm.Greens, len(species_list), 2, 1)
	for i in range(len(species_list)):
		species_color_map[species_list[i]] = m.to_rgba(i)

	# Special: set '-' to gray
	species_color_map['-'] = (0.6, 0.6, 0.6, 1.0)

	return species_color_map, ordered_species_list

def get_cm_ScalerMappable(cmap, num_colors, offset1 = 0, offset2 = 0):
	norm = matplotlib.colors.Normalize(vmin = 0 - offset1, vmax = num_colors- 1 + offset2)
	return matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

def list_to_colors(input_list):
	cmap = matplotlib.cm.hsv

	input_list_dict = {}
	i = 0
	for elem in set(input_list):
		input_list_dict[elem] = i
		i += 1

	norm = matplotlib.colors.Normalize(vmin=0, vmax=i-1)
	m = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

	color_list = []
	for elem in input_list:
		color_i = input_list_dict[elem]
		color_list.append(m.to_rgba(color_i))

	return color_list

def colors_to_legend_elements(colors, labels):
	legend_elements = {}
	for color, label in zip(colors, labels):
		legend_elements[color] = matplotlib.patches.Patch(facecolor=color, label=label)
	return [legend_elements[c] for c in remove_list_duplicates(colors)]

def remove_list_duplicates(mylist):
	# Returns list with only unique elements (like set)
	# but ordered according to first appearnce in original list
	result = []
	for elem in mylist:
		if elem not in result:
			result.append(elem)
	return result
