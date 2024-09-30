def add_gene_programs(adata):
  """ After running cNMF add each gene program to the anndata object
  Args:
  adata: the anndata object that cNMF was run on
  """
    program_names = []
    # go through all the program names and make an array
    for value in usage_norm.columns:
        program_names.append(str(value))
    usage_norm.columns = program_names
    # create a new obs for each gene program and add the gene program expression
    data.obs[usage_norm.columns] = usage_norm
