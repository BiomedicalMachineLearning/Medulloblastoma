import numpy as np
from anndata import AnnData

def species_split(data, species='human', remove_genes=False,
                  out_format=None):
    """ Splits anndata into data/mix spots for data genes.
        remove_genes indicates whether to remove mitochondrial/ribosomal genes.
        out_format='human' or 'mouse' to format genes either way, default is same
        as input species.
    """
    species_labels = data.obs['species'].values.astype(str)
    mix_bool = species_labels == 'mix'
    human_bool = species_labels == 'human'
    mouse_bool = species_labels == 'mouse'

    h_prefix, m_prefix = 'hg38-', 'mm10-'
    prefix = h_prefix if 'human' in species else m_prefix
    if species == 'human':
        species_bool = np.logical_or(human_bool, mix_bool)
    elif species == 'mouse':
        species_bool = np.logical_or(mouse_bool, mix_bool)
    elif species == 'human-only':
        species_bool = human_bool
    elif species == 'mouse-only':
        species_bool = mouse_bool
    elif species == 'mix-only':
        species_bool = mix_bool
        prefix = ''
    elif species == 'mix-humang':
        species_bool = mix_bool
    elif species == 'mix-mouseg':
        species_bool = mix_bool

    # NOT implimented.... could be.
    if remove_genes and 'human' in species:
        remove_prefixes = ['MT-', 'RPL', 'RPS']
    elif remove_genes and 'mouse' in species:
        remove_prefixes = ['mt-', 'Rpl', 'Rps']
    elif remove_genes and 'mix-only':
        remove_prefixes = ['MT-', 'RPL', 'RPS', 'mt-', 'Rpl', 'Rps']
    else:
        remove_prefixes = []

    data_species = data[species_bool, :]

    species_genes = [gene for gene in data.var_names if prefix in gene]
    species_df = data_species[:,species_genes].to_df()
    split_cond = prefix if prefix!='' else ' '
    species_genes_formatted = np.array([gene.split(split_cond)[-1]
                                        for gene in species_genes])
    species_df.columns = species_genes_formatted

    if out_format=='human' and prefix!='hg38-':
        genes = [gene.upper() for gene in species_df.columns]
        species_df.columns = genes
    elif out_format=='mouse' and prefix!='mm10-':
        genes = [gene[0].upper()+gene[1:].lower() for gene in species_df.columns]
        species_df.columns = genes

    species_ad = AnnData(species_df, obs=data_species.obs, uns=data_species.uns,
                                 obsm=data_species.obsm, obsp=data_species.obsp)

    return species_ad







