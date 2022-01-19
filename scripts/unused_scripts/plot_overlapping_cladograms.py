from __future__ import division
import os
import ete3

from itertools import combinations

import toytree
import toyplot
import numpy


directory = '/Users/williamrshoemaker/GitHub/snv_prevalence/'

tree_directory = '%sdata/cladogram_newick/' % directory


mtree_list = []

# get overlapping hosts

host_count_dict = {}

for gene_file in os.listdir(tree_directory):
    if gene_file.endswith(".txt") == False:
        continue
    gene_path = '%s%s' % (tree_directory, gene_file)
    gene_tree = ete3.Tree(gene_path, 1)

    species_name = gene_file.split('.')[0]
    node_labels = [t_i.name for t_i in gene_tree.get_leaves()]

    for node_label in node_labels:
        if node_label not in host_count_dict:
            host_count_dict[node_label] = []

        host_count_dict[node_label].append(species_name)





host_count_dict_to_keep = {k: v for k, v in host_count_dict.items() if len(v) >= 10}

hosts_to_keep = list(host_count_dict_to_keep.keys())


hosts_to_keep_pairs = combinations(hosts_to_keep, 2)

print(len(hosts_to_keep))

jaccar_list = []
intersection_list = []
for pair in hosts_to_keep_pairs:

    host_i_set = set(host_count_dict_to_keep[pair[0]])
    host_j_set = set(host_count_dict_to_keep[pair[1]])

    jaccard_pair = len(host_i_set & host_j_set) / len(host_i_set | host_j_set)
    intersection_list.append(len(host_i_set & host_j_set))

    jaccar_list.append(jaccard_pair)

print(numpy.mean(jaccar_list))

print(numpy.mean(intersection_list))
#fix this later

for gene_file in os.listdir(tree_directory):
    if gene_file.endswith(".txt") == False:
        continue
    gene_path = '%s%s' % (tree_directory, gene_file)
    gene_tree = ete3.Tree(gene_path, 1)

    ec = gene_file.split('.')[0]

    # check for polytomies
    polytomy_status = False
    for n in gene_tree.get_descendants():
        if len(n.children) > 2:
            polytomy_status = True

    if polytomy_status == True:
        continue


    tre = toytree.tree(gene_path)

    #print(tre.mod.make_ultrametric())

    tre = tre.mod.make_ultrametric()

    mtree_list.append(tre)




mtree = toytree.mtree(mtree_list)


for idx, tree in enumerate(mtree.treelist):

    tree.style.edge_style["stroke"] = toyplot.color.brewer.palette("BlueGreen")[4]
    tree.style.edge_style['stroke-opacity'] = 0.05


#fixed_order=species_tree.get_tip_labels(),

# draw a cloud tree which enforces a fixed tip order
canvas, axes, mark = mtree.draw_cloud_tree()
#canvas, axes, mark = mtree.draw_cloud_tree(
#    tip_labels_style={"font-size": "11px"},
#    height=750,
#    width=450,
#);




# write as PDF
import toyplot.pdf
toyplot.pdf.render(canvas, "%sanalysis/cloudtree.pdf" % directory)
