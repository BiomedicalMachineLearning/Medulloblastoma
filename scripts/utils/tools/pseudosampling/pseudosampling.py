""" Goal of this script is to impliment useful method for pseudobulking,
    such as random-sampling within a particular condition & taking mean, sum,
    or median. Useful for doing DE in a way that alleviates the pseudoreplication
    problem.
"""

import numpy as np
import pandas as pd

import scripts.utils.utils.inputHelpers as ihs

def pseudo_sample(expr, sample_labels, k, method='mean', sample_set=None):
    """Performs pseudosampling, whereby within defined groups random samples \
        are taken, then averaged, to get the final estimate.

    Args:
        expr (pd.DataFrame): Genes*Cells data frame.
        sample_labels (np.array): The labels within which data will be pseudo \
                             bulked by subsampling random indices & aggregating.
        sample_set (np.array): The unique labels of the set, if input can be \
                            useful to control output order, otherwise will \
                            just use np.unique(sample_labels).
        k (int): The number of subsamples to take.
        method (str): The technique to use for aggregation; either \
                                                           mean, median, or sum.

    Returns:
         pd.DataFrame: Genes*Pseudosamples, where number of pseudosamples is \
                        k*len(sample_set).
    """
    sample_set = sample_set \
                   if type(sample_set)!=type(None) else np.unique(sample_labels)

    metaspot_scores = np.zeros((expr.shape[0], k*len(sample_set)))
    sample_i = 0
    n_spots = []
    pseudolabels = []
    for i, sample in enumerate(sample_set):
        sample_indices = np.where(sample_labels==sample)[0]
        n_rand_size = len(sample_indices)//k
        for kth_sample in range(k-1):
            rand_indices = np.random.choice(sample_indices, n_rand_size,
                                            replace=False)
            n_spots.append( len(rand_indices) )
            sample_indices = [index for index in sample_indices
                                     if index not in rand_indices]
            metaspot_scores[:, sample_i] = get_bulked(
                                           expr.values[:,rand_indices], method)
            pseudolabels.append(f'{sample}_{kth_sample}')
            sample_i += 1

        metaspot_scores[:, sample_i] = get_bulked(expr.values[:,sample_indices],
                                                                         method)
        n_spots.append( len(sample_indices) )
        pseudolabels.append(f'{sample}_{k}')
        sample_i += 1

    metaspot_scores = pd.DataFrame(metaspot_scores,
                                   index=expr.index, columns=pseudolabels)

    return metaspot_scores

def get_bulked(values, method):
    """Performs the summarising of the values, either by mean, median, or sum."""
    ihs.raiseErrorIfNotPresent(method, ['mean', 'media', 'sum'])

    if method == 'mean':
        return values.mean(axis=1)
    elif method == 'median':
        return values.median(axis=1)
    else:
        return values.sum(axis=1)








