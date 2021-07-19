# BrainCancer
For sharing scripts/data on analyses for the MedulloBlastoma data.

## Quick Index

    docs/ -> Contains project documentation.
        index.md -> Contains the index which describes each script in terms of 
                    INPUT/OUTPUT, and brief description of function.
                    
    scripts/ -> Contains the scripts described in docs/index.md.

    data/
        ryan_ids/ -> Contains the spot barcodes that were left after very 
                        basic/lenient filtering (>200genes)
                        
# Figure 4 (19th July 2021)
The figure 4 scripts can be found in:

    scripts/species_classify_v2/mix_only/
    
Script descriptions:

        mix_analysis.py (07/07/2021) -> Purpose of this script is to load only the
                         mix spots, but separate into human genes/mouse genes.
                         For the mix-humang, compare to Hovestadt reference.
                         For the mix-mouseg, compare to Vladoiu reference. 
                         Going to do this for each sample independently, but then
                         add the scores from the comparison to the origin mix
                         anndata, & subsequently integrate the mix spots into a 
                         common space. Thereafter, will examine the scores &
                         cluster, labelling the clusters according to the 
                         cell type combinations. 
                         
                         INPUT: data/Visium8*/ 
                                data/thirdy_party_data/Hovestadt2019_Nature_scRNA/
                                                                  hov_scrna.h5ad
                                data/thirdy_party_data/Vladoiu2019_Nature_scRNA/
                                                        Vladoiu2019_logcpms.h5ad
 
                         OUTPUT: figure_components/species_classify_v2/mix_only/
                                 data/scanpy_h5ads/integrated_border.h5ad
                                 
                                 data/thirdy_party_data/Vladoiu2019_Nature_scRNA/
                                                        Vladoiu2019_lateRef.h5ad
                                                        
    border_figures.py (08/07/2021) -> 
                        Based on the analysis conducted in mix_analysis.py, 
                        loads in the data & generates the figure panels to explain
                        the observations. 
                        
                        INPUT: data/scanpy_h5ads/integrated_border.h5ad
                               data/thirdy_party_data/Hovestadt2019_Nature_scRNA/
                                                                  hov_scrna.h5ad
                                data/thirdy_party_data/Vladoiu2019_Nature_scRNA/
                                                        Vladoiu2019_lateRef.h5ad
                                                        
                        OUTPUT: figure_components/species_classify_v2/
                                                               mix_only/panels/*     




