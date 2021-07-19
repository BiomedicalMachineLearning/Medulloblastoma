""" Created: 19th April 2021

Using the stlearn clusters, normalises/scales for each spot & then performs the
rank-ordering followed by FET enrichment for each gene above a particular cutoff
within each cluster of spots
(comparing frequency about cutoff in list versus below cutoff).
perhaps impliment automatic detection of a good cutoff based on inflection point)

OUTPUT: data/DE_out/stlearn_cluster_DE/triage_spot_profiles.xlsx
"""


################################################################################
                    # Environment setup #
################################################################################

import os
import stlearn as st
import scanpy as sc
import pandas as pd
import pandas
import numpy

import scripts.helpers as hs
hs.setUp()
import scripts.stlearn_helpers as st_hs

import beautifulcells.tools.triage.triage_analysis as triage_analysis
import beautifulcells.postprocessing.format.format as format

################################################################################
                    # Loading data & gene filtering #
################################################################################
data_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/data/'
samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
data_dirs = [data_dir+sample for sample in samples_]
species_metas = ['A1_treated_species.txt', 'B1_treated_species.txt',
                     'C1_untreated_species.txt', 'D1_untreated_species.txt']

filt_files = numpy.array( os.listdir(data_dir+'ryan_ids/') )
filt_files = filt_files[numpy.argsort(filt_files)]
spot_filters = [pd.read_csv(data_dir+'ryan_ids/'+filt_file,
                           header=None).values[:,0] for filt_file in filt_files]

samples = ['A1', 'B1', 'C1', 'D1']

# Loading in the un-normalised counts, filtered to data/mix spots.
datas = [st_hs.load(sample_dir, spot_filters[i], species_metas[i])
         for i, sample_dir in enumerate(data_dirs)]
# Removing lowly expressed genes
[st.pp.filter_genes(data, min_cells=5) for data in datas]

################################################################################
                    # Pseudobulk & normalise clusters #
################################################################################
# Loading the cluster labels #
meta = pandas.read_csv(data_dir+'spot_meta/clusters_human_mix.txt',
                                                          sep='\t', index_col=0)
tissue = meta.loc[:,'tissue_type'].values

sheet_names = []
triages = []
triage_objs = []
for i, data in enumerate(datas):
    sample = samples[i]
    sample_indices = [i for i in range(meta.shape[0])
                      if sample in meta.index.values[i]]
    ordered_obs = meta.index.values[sample_indices]
    obs_names = numpy.array([f'{sample}-{name}' for name in data.obs_names])
    counts = data.to_df().transpose()
    counts.columns = obs_names
    counts = counts.loc[:,ordered_obs]
    print(numpy.all(counts.columns.values==ordered_obs))

    # Subsetting to data genes #
    human_genes = hs.getSpeciesGenes(counts, 'hg38-')
    counts = counts.loc[human_genes, :]
    counts.index = hs.getSpeciesGenes(counts, 'hg38-', remove_prefix=True)

    # Normalising & Scaling #
    normed = st_hs.norm_df(counts, scale=True)

    # Performing TRIAGE analysis #
    sample_tissues = tissue[sample_indices]
    triage = triage_analysis.Triage(organism='data', average=False,
                                                                 binary_data='')
    triaged = triage.fit_transform(counts, labels=sample_tissues)
    triage_profiles = triage.test_enrichment()

    sheet_names.append(f'{sample}_TRIAGED')
    triages.append( triage_profiles )
    triage_objs.append( triage )

format.writeDFsToExcelSheets(data_dir+
                          'DE_out/stlearn_cluster_DE/triage_spot_profiles.xlsx',
                                                           triages, sheet_names)

"""
# Testing idea of using count distributions #
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = triage_objs[3]._label_conts['proliferative-stem-like']
props = data.values[0,:]/(data.values[0,:]+data.values[1,:])
print(len(props))
# 1046 for 0 pairs,
# plt.hist(data.values[0,:], bins=100, density=True)
# plt.show()
#
# plt.hist(data.values[0,:]/data.shape[1], bins=100, density=True)
# plt.show()

# Cutoff approach #
order = numpy.argsort(-props)
props_ordered = props[order]

plt.scatter(list(range(len(props_ordered))), props_ordered)
plt.show()

total = sum(props_ordered)
dens = [props_ordered[i]/total for i in range(1, len(props_ordered))]

plt.scatter(list(range(len(dens))), dens)
plt.show()

diffs = numpy.diff(props_ordered)
#diffs = numpy.diff(diffs)

# Fitting the interpolation #
from scipy.interpolate import interp1d

f = interp1d(list(range(len(diffs))), diffs, 'linear')
diff_infer = [f(i/3) for i in range(len(diffs))]

n = 1000
plt.scatter(list(range(len(diff_infer)))[0:n], diff_infer[0:n])
plt.show()

n = 1000
plt.scatter(list(range(len(diffs)))[0:n], diffs[0:n])
plt.show()

def f(x,
	  a, b, k, x0, c, #decay
	):
    decay = (b * (1 + numpy.exp(-k * (x)))) +x0 #+ (c*x) + x0
    #decay = a + (b*numpy.exp(-k*x)) + (c*x)
    #decay = numpy.array(decay)
    #decay[decay < 0] = 0
    return decay

opt_parms, parm_cov = curve_fit(f,
                                list(range(len(diffs)))[0:n],
                                diff_infer[0:n],
                                #p0=(0, -.1, 0.5, 0.6),
                                maxfev=100000)

sim = numpy.array( [f(val, *opt_parms) for val in range(n)] )

plt.scatter(list(range(n)), sim)
plt.show()

bool_ = numpy.round(sim, 4) == numpy.round(opt_parms[1]+opt_parms[3], 4)
cutoff = numpy.where(bool_)[0][0]
print(cutoff)

# Trying Pareto distribution #
from scipy.stats import pareto

total = sum(props_ordered)
dens = [props_ordered[i]/total for i in range(1, len(props_ordered))]

plt.scatter(list(range(len(dens))), dens)
plt.show()

params = pareto.fit(props_ordered)
vals = pareto.rvs(params[0], loc=params[1], scale=params[2], size=1000)

order = numpy.argsort(-vals)
vals = vals[order]

plt.scatter(list(range(len(vals))), vals)
plt.show()

plt.scatter(list(range(len(vals))), props_ordered[0:len(vals)])
plt.show()

pareto.cdf(5, params[0], loc=params[1], scale=params[2])

#### Pareto distribution fit #####
def f(x,
	  a, b, k, x0, c, #decay
	):
    decay = b / ((x-x0)**(b+1))
    decay = decay / c
    return decay

opt_parms, parm_cov = curve_fit(f,
                                list(range(len(diffs)))[0:n],
                                diff_infer[0:n],
                                #p0=(0, -.1, 0.5, 0.6),
                                maxfev=100000)

sim = numpy.array( [f(val, *opt_parms) for val in range(n)] )

plt.scatter(list(range(n)), sim)
plt.show()


# negative binomial #s
import scipy
import numpy as np
import statsmodels.api as sm

background = data.values[0,:]
scores = background

pmin, pmax = min(background), max(background)
background2 = [item - pmin for item in background]
x = np.linspace(pmin, pmax, 1000)
res = sm.NegativeBinomial(
    background2, np.ones(len(background2)), loglike_method="nb2"
).fit()#start_params=[0.1, 0.3])
mu = np.exp(res.params[0])
alpha = res.params[1]
Q = 0
size = 1.0 / alpha * mu ** Q
prob = size / (size + mu)

from scipy.stats import nbinom
from statsmodels.stats.multitest import multipletests

r = nbinom.rvs(size, prob, size=data.shape[1])
plt.hist(r, bins=100, density=True)
plt.show()

plt.hist(data.values[0,:], bins=100, density=True)
plt.show()

# Calculate probability for all spots
pvals = 1 - scipy.stats.nbinom.cdf(scores - pmin, size, prob)
pval_adj = multipletests(pvals, 0.05, 'fdr_bh')[1]

print( len(numpy.where((pvals)<0.01)[0]) )
print( len(numpy.where((pval_adj)<0.05)[0]) )

print(data.columns.values[pvals<.01])

# Poisson #
pmin, pmax = min(background), max(background)
background2 = [item - pmin for item in background]
x = np.linspace(pmin, pmax, 1000)
res = sm.Poisson(background2, np.ones_like(background2)).fit()

# beta distribution #
from scipy.stats import beta

props = background/len(background)
a1, b1, loc1, scale1 = beta.fit(props)

r2 = beta.rvs(a1, b1, size=data.shape[1])
plt.hist(r2, bins=100, density=True)
plt.show()

pvals = 1 - scipy.stats.beta.cdf(props, a1, b1)
pval_adj = multipletests(pvals, 0.05, 'bonferroni')[1]

print( len(numpy.where((pvals)<0.05)[0]) )
print( len(numpy.where((pval_adj)<0.05)[0]) )
"""

