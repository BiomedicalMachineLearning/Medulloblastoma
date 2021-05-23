""" Created: March 24th 2021
    Conda env: cell2location

Ref:
https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo_downstream.html

The purpose of this script is to take as input a dataframe of cell proportion 
estimates across spots, and subsequently perform NNMF using cell2location on these 
to identify 'tissue types' - or factors that represent combinations of cells.
"""

import os, sys
sys.path.append('/home/uqbbalde/myPython/BeautifulCells/')
import beautifulcells.visualisation.helpers as vhs
import numpy, pandas

data_type = 'float32'
import cell2location.models as c2l

import matplotlib as mpl
from matplotlib import pyplot as plt

import warnings
warnings.filterwarnings('ignore')

################################################################################
                        # Loading data #
################################################################################
base_dir = '/30days/uqbbalde/MedullaBlastoma/'
plot_dir = base_dir+'plots/decon_out/'
decon_dir = base_dir+'data/decon_out/'
decon_file = 'human_vladoiu_cellprops.txt'

decon_df = pandas.read_csv(decon_dir+decon_file, sep='\t', 
			   index_col=0, header=0)
decon_df = decon_df.iloc[0:decon_df.shape[0]-1,:]

zero_bool = decon_df.values == 0
cell_bool = numpy.apply_along_axis(numpy.all, 1, zero_bool)==False

decon_df = decon_df.loc[cell_bool,:]

decon_array = decon_df.values.transpose()

################################################################################
                        # Performing the NNMF #
################################################################################
n_fact = int(10)

# create model class
mod_sk = c2l.CoLocatedGroupsSklearnNMF(n_fact, decon_array,
        n_iter = 10000,
        verbose = True,
        var_names=decon_df.index.values,
        obs_names=decon_df.columns.values,
        fact_names=['fact_' + str(i) for i in range(n_fact)],
        sample_id=None,
        init='random', random_state=0,
        nmf_kwd_args={'tol':0.0001})

# train 5 times to evaluate stability
mod_sk.fit(n=5, n_type='restart')

########## Model checking plots ############
# Looking at the concordance between the factors using different
# optimisation starts.
with mpl.rc_context({'figure.figsize': (10, 8)}):
    mod_sk.evaluate_stability('cell_type_factors', align=True)

vhs.dealWithPlot(True, False, True, plot_dir, 
		'labelTransfer_human_vladoiu_cell2LocFactorConcord.pdf', 300)

# evaluate accuracy of the model
mod_sk.compute_expected()
mod_sk.plot_posterior_mu_vs_data()

vhs.dealWithPlot(True, False, True, plot_dir,
                'labelTransfer_human_vladoiu_cell2LocFactorImpute.pdf', 300)

########## Visualising results ############
mod_sk.sample2df(node_name='nUMI_factors', ct_node_name = 'cell_type_factors')

# print the fraction of cells of each type located to each combination
top_loads = mod_sk.print_gene_loadings(loadings_attr='cell_type_fractions',
           		              gene_fact_name='cell_type_fractions')
print(top_loads)

top_loads.to_csv(decon_dir+'labelTransfer_human_vladoiu_cell2LocFactorTopLoadings.txt',
		sep='\t')

# Plotting a heatmap of the loadings #
from re import sub
mod_sk.cell_type_fractions.columns = [sub('mean_cell_type_factors', '', i)
                                      for i in mod_sk.cell_type_fractions.columns]

mod_sk.plot_gene_loadings(mod_sk.var_names_read, mod_sk.var_names_read,
                        fact_filt=mod_sk.fact_filt,
                        loadings_attr='cell_type_fractions',
                        gene_fact_name='cell_type_fractions',
                        cmap='RdPu', figsize=[10, 15])

vhs.dealWithPlot(True, False, True, plot_dir,
                'labelTransfer_human_vladoiu_cell2LocFactorLoadings.pdf', 300)

###### Saving the model ######
# save co-location models object
def pickle_model(mod, path, file_suffix=''):
    file = path + 'Model_' + str(mod.__class__.__name__) + '_' + str(mod.n_fact) + '_' + file_suffix
    pickle.dump({'mod': mod, 'fact_names': mod.fact_names}, file = open(file, "wb"))
    print(file)

pickle_model(mod_sk, decon_dir+'labelTransfer_human_vladoiu_cell2Loc', file_suffix='.pkl')

###### Saving the factor scores for each spot #########
mod_sk.location_factors_df.to_csv(decon_dir+'labelTransfer_human_vladoiu_cell2LocScores.txt',
				sep='\t')



