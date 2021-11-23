
def get_pretty_species_name(species_name, include_number=False):

    items = species_name.split("_")

    pretty_name = "%s %s" % (items[0], items[1])

    if include_number:
        pretty_name += (" (%s)" % (items[2]))

    return pretty_name

def get_abbreviated_species_name(species_name):

    items = species_name.split("_")

    pretty_name = "%s. %s" % (items[0][0], items[1])

    return pretty_name


clade_type_label_dict = {'all': 'All hosts', 'largest_clade': 'Largest clade'}


sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i']
