import os
import pickle



for filename in os.listdir('/Users/wrshoemaker/GitHub/negative_selection_microbiome/data/strain_data_3/'):
    if filename.endswith(".pkl"):
        # print(os.path.join(directory, filename))
        old_path = '/Users/wrshoemaker/GitHub/negative_selection_microbiome/data/strain_data_3/' + filename
        with open(old_path, 'rb') as handle:
            b = pickle.load(handle)

        output_name = '/Users/wrshoemaker/GitHub/negative_selection_microbiome/data/strain_data/' + filename
        with open(output_name, 'wb') as handle:
            pickle.dump(b, handle,  protocol=2)

    else:
        continue
