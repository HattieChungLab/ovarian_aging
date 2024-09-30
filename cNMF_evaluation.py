# After running cNMF this file can be run to inspect the number of gene programs that cNMF recommends
# and the top genes in each program 

import pickle
from IPython.display import Image, display

# name of the file cNMF was ran on
path = 'your_file_name'
pathToPickle = path + '.cnmf_run_orig_counts.pkl'
with open(pathToPickle, 'rb') as f:
    cnmf_obj = pickle.load(f)

# will show a graph with stability and error at each value of k, typically the k value
# is chosen by the highest stability as error will always increase as k increases, but
# this can lead to overfitting 
display(Image('cNMF_outputs/ + path + '/' + path + '.k_selection.png'))

# selected_K value based on visual evaluation from (cNMF_file_name).k_selection.png
selected_K = 8
density_thresh = 0.1
cnmf_obj.consensus(k=selected_K, density_threshold=density_thresh)

usage_norm, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(K=selected_K, 
                                                                           density_threshold=density_thresh)

# save a csv file of the top 10 genes in each program for further analyses
top_genes_10 = top_genes.head(10)
top_genes_10.to_csv('cNMF_top_10_genes.csv')
