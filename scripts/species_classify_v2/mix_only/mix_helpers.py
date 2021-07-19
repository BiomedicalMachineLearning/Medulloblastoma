"""
    Helper functions for the mix only analysis.
"""

import numpy as np

def vlaodiu_label_merge(cell_labels, label_set):
    """ Merges the vladoiu labels for simpler label transfer !
    """
    merged = {'glial': ['Astrocyte/Bergmann glia precursors',
                        'Oligodendrocyte precursor cells',
                       'Gliogenic progenitors-1', 'Gliogenic progenitors-2'],
             'Granule cells': ['Embryonic and postnatal GCPs-1',
                               'Granule cells',
                               'Postnatal GCPs-2'],
             'vasculature': ['Endothelial cells', 'Meninges', 'Pericytes'],
             'GABA interneurons': ['GABA interneurons',
                                   'GABA interneuron precursors'],
             'Microglia': ['Microglia'],
             'Glutamatergic neurons': [
                 'Postnatal excitatory cerebellar nuclei neurons'],
             'Unipolar brush cells': ['Unipolar brush cell precursors',
                                      'Unipolar brush cells'],
             'Purkinje cells': ['Purkinje cells']
             }

    # Checking to make sure everything named correctly & accounted for #
    accounted = []
    for group in merged:
        for cell_type in merged[group]:
            accounted.append(cell_type)
            if cell_type not in label_set:
                print(cell_type)
    print('Missing cell types:')
    print([cell_type for cell_type in label_set if
           cell_type not in accounted])

    #Inversing the mapping #
    rev_merge = {}
    for group in merged:
        for cell_type in merged[group]:
            rev_merge[cell_type] = group

    new_labels = np.array([rev_merge[cell_label] for cell_label in cell_labels])

    return new_labels
























