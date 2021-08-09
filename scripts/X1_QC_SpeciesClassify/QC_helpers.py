""" QC helper functions
"""

def getGenePercents(data, genes, n_expr=True):
    """n_expr=True if in terms of boolean numbers, or in-terms of total counts \
    per spot to calculate the percentage.
    """

    df = data.to_df()
    df_sub = df.loc[:, genes]

    if n_expr:
        n_genes = (df.values > 0).sum(axis=1)
        n_sub = (df_sub.values > 0).sum(axis=1)
        perc = n_sub/n_genes

    else:
        n_genes = (df.values).sum(axis=1)
        n_sub = (df_sub.values).sum(axis=1)
        perc = n_sub / n_genes

    return perc


