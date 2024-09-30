import sys
#cNMF_dir = 'adataframes'
#name = 'orig_counts_adata_raw_2groups_pooled_young_endothelial'

cNMF_dir = sys.argv[1]
name = sys.argv[2]

def run_cNMF(data_path, run_name, numhvgenes=2000):
    from cnmf import cNMF
    import os

    #cNMF_dir = 'ovarian_aging/analysis/segment_analyses/follicles'
    cNMF_outdir = 'cNMF_outputs'
    
    numiter=100 # Number of NMF replicates. Set this to a larger value ~200 for real data. We set this to a relatively low value here for illustration at a faster speed
    numhvgenes=numhvgenes ## Number of over-dispersed genes to use for running the actual factorizations
    K_range = range(5,25)

    ## Results will be saved to [output_directory]/[run_name] which in this example is example_PBMC/cNMF/pbmc_cNMF
    output_directory = '%s' %cNMF_outdir
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    ## Specify the Ks to use as a space separated list in this case "5 6 7 8 9 10"
    K = ' '.join([str(i) for i in K_range])

    seed = 14 ## Specify a seed pseudorandom number generation for reproducibility

    ## Path to the filtered counts dataset we output previously

    ## Initialize the cnmf object that will be used to run analyses
    cnmf_obj = cNMF(output_dir=output_directory, name=run_name)

    ## Prepare the data, I.e. subset to 2000 high-variance genes, and variance normalize
    cnmf_obj.prepare(counts_fn=data_path, components=K_range, n_iter=numiter, seed=14, num_highvar_genes=numhvgenes)

    ## Specify that the jobs are being distributed over a single worker (total_workers=1) and then launch that worker
    cnmf_obj.factorize(worker_i=0, total_workers=1)

    cnmf_obj.combine()
    cnmf_obj.k_selection_plot(close_fig=False)
    print('This saves the corresponding figure to the following file: %s' % cnmf_obj.paths['k_selection_plot'])
    
    return cnmf_obj

def get_corr_df(ad, module_genes): 
    if isinstance(module_genes, dict): 
        top_module_genes = list(set([item for sublist in module_genes.values() for item in sublist]))
    else: 
        top_module_genes = module_genes
    DEG_df = pd.DataFrame(ad[:,top_module_genes].layers['processed'], 
                          index=ad.obs.index, columns=top_module_genes)
    DEG_corr = DEG_df.corr()
    return DEG_corr

def calculate_pvalues(df):
    from scipy.stats import pearsonr
    df = df.dropna()._get_numeric_data()
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            pvalues[r][c] = round(pearsonr(df[r], df[c])[1], 4)
    return pvalues

def get_mask_subbed_df(df, mask_val=-100): 
    # substitute with mask_val 
    mask = df.isnull()
    subbed_df = df.fillna(mask_val)
    return subbed_df

import seaborn as sns

def genes_df_to_list(gene_list, n_top=20):
    mod_genes = gene_list.head(n_top).values
    top_module_genes = list(set([item for sublist in mod_genes for item in sublist]))
    return top_module_genes

def make_corr_df(gene_list, ad, keep_sig=False, n_top=20): 
    
    corr_thresh=0.05
    mask_val = 0
    top_module_genes = genes_df_to_list(gene_list, n_top=n_top)
    
    DEG_df = pd.DataFrame(ad[:,top_module_genes].layers['processed'], 
                      index=ad.obs.index, columns=top_module_genes)
    DEG_corr = DEG_df.corr()
    DEG_corr.fillna(0, inplace=True)
    
    if keep_sig: 
        df_pvals = calculate_pvalues(DEG_corr)
        df_corr_sig = DEG_corr[df_pvals<corr_thresh]
        df_corr_masked = get_mask_subbed_df(df_corr_sig, mask_val)

        num_ele = df_corr_masked.shape[0]
        genes_to_remove_0 = ((df_corr_masked==0).sum(axis=0)==(num_ele-1))
        genes_to_remove_1 = ((df_corr_masked==0).sum(axis=1)==(num_ele-1))

        remove_genes = genes_to_remove_0[genes_to_remove_0 & genes_to_remove_1].index

        df_corr_masked.drop(remove_genes, axis=0, inplace=True)
        df_corr_masked.drop(remove_genes, axis=1, inplace=True)
        df_plot = df_corr_masked
        
    else:
        df_plot = DEG_corr
        
    sns.set(font_scale=0.6)
    g = sns.clustermap(df_plot,
                    xticklabels=[], 
                    yticklabels=DEG_corr.index,
    #                 col_colors=genes_color,
                    vmin=-0.75, vmax=0.75, 
                    cmap='bwr',
                    dendrogram_ratio=0.08, 
                    colors_ratio=0.025,
                    cbar_pos=[1,0.6,0.02,0.1])
    g.ax_col_dendrogram.set_visible(False)
    g.ax_row_dendrogram.set_visible(False)
    gene_order = g.dendrogram_row.reordered_ind
    ax = g.ax_heatmap
    ax.axhline(y=0, color='k',linewidth=2)
    ax.axvline(x=0, color='k',linewidth=2)
    sns.reset_orig()
    return [top_module_genes[g] for g in gene_order]

def make_ordered_corr(ordered_gene_list, ad, keep_sig=False): 
    
    plt.figure()
    
    corr_thresh=0.05
    mask_val = 0
    
    DEG_df = pd.DataFrame(ad[:,ordered_gene_list].layers['processed'], 
                      index=ad.obs.index, columns=ordered_gene_list)
    DEG_corr = DEG_df.corr()
    DEG_corr.fillna(0, inplace=True)
    
    if keep_sig: 
        df_pvals = calculate_pvalues(DEG_corr)
        df_corr_sig = DEG_corr[df_pvals<corr_thresh]
        df_corr_masked = get_mask_subbed_df(df_corr_sig, mask_val)

        num_ele = df_corr_masked.shape[0]
        genes_to_remove_0 = ((df_corr_masked==0).sum(axis=0)==(num_ele-1))
        genes_to_remove_1 = ((df_corr_masked==0).sum(axis=1)==(num_ele-1))

        remove_genes = genes_to_remove_0[genes_to_remove_0 & genes_to_remove_1].index

        df_corr_masked.drop(remove_genes, axis=0, inplace=True)
        df_corr_masked.drop(remove_genes, axis=1, inplace=True)
        df_plot = df_corr_masked
        
    else:
        df_plot = DEG_corr
        
    sns.set(font_scale=0.6)
    g = sns.heatmap(df_plot,
                    xticklabels=[], 
                    yticklabels=DEG_corr.index,
    #                 col_colors=genes_color,
                    vmin=-0.75, vmax=0.75, 
                    cmap='bwr',
                   square=True,
                   cbar=False)
#     g.ax_col_dendrogram.set_visible(False)
#     g.ax_row_dendrogram.set_visible(False)
#     gene_order = g.dendrogram_row.reordered_ind
#     ax = g.ax_heatmap
#     ax.axhline(y=0, color='k',linewidth=2)
#     ax.axvline(x=0, color='k',linewidth=2)
    sns.reset_orig()



cnmf_obj = run_cNMF('/home/ss4848/{}/{}'.format(cNMF_dir, name), name)






import pickle

with open('{}.cnmf_run_orig_counts.pkl'.format(sys.argv[2]), 'wb') as f:
    pickle.dump(cnmf_obj, f)
