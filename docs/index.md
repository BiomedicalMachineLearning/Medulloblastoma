# Index
Index which describes each script in scripts/ in terms of INPUT/OUTPUT, and 
brief description of function.

X1_QC_SpeciesClassify/

    X1_QC.py -> Performs the QC to filter spots & looks at the distribution of 
                                                               the counts/genes.
             INPUT:  * data/Visium8_{sample_name}_Hybrid/*
             OUTPUT: * data/filter_ids/{sample_name}_filtered.txt -> 
                              Contains the barcodes of the spots that passed QC.       
                     * figure_components/QC_figure/*
                     
    X2_generate_seurat_rds.R -> Generates seurat rds with sctransform normalised
                                data which is later used for species scoring &
                                giotto enrichment analysis.
                                NOTE: memory intensive, had to run on a hpc.
                                
                    INPUT: * data/Visium8_{sample_name}_Hybrid/*
                           * data/filter_ids/{sample_name}_filtered.txt
                           
                    OUTPUT: * seurat_rds/*
                    
    X3_human_mouse_score.R -> Sums to sctransformed normalised gene expression
                             for human/mouse genes to generate human/mouse 
                             scores per spot in each sample, to subsequently
                             allow for species clasification.
                             
                    INPUT: * seurat_rds/*treated.rds
                    OUTPUT: * spot_meta/*_species.txt
                    
    X4_species_classify.py -> Uses the species scores to classify spots by into
                                                                human/mix/mouse.
                                    
                     INPUT: * data/Visium8*/
                            * spot_meta/*_species.txt
                     OUTPUT: * data/spot_meta/species_classify_v2/
                             * figure_components/species_figure/
                             
X2_DEAnalysis/

    X1_Pseudo_Limma_DE_human.R -> Performs the DE between Palbociclib treated & 
                                 control samples pseudobulked human spots/genes.
                            
                            INPUT: * data/seurat_rds/all.rds
                                   * data/spot_meta/species_classify_v2/
                                                                   *_species.txt
                            
                            OUTPUT: * data/DE_out/Pseudo_Limma_human/*
                            
    X2_Pseudo_Limma_DE_mixHuman.R -> Performs the DE between Palbociclib treated 
                           & control samples pseudobulked mix spots/human genes.
    
                            INPUT: * data/seurat_rds/all.rds
                                   * data/spot_meta/species_classify_v2/
                                                                   *_species.txt
                            
                            OUTPUT: * data/DE_out/Pseudo_Limma_mixHuman/*
    
    X3_Pseudo_Limma_DE_mixMouse.R -> Performs the DE between Palbociclib treated 
                           & control samples pseudobulked mix spots/mouse genes.
                           
                           INPUT: * data/seurat_rds/all.rds
                                  * data/spot_meta/species_classify_v2/
                                                                   *_species.txt
                            
                           OUTPUT: * data/DE_out/Pseudo_Limma_mixMouse/*
    
    X4_Pseudo_Limma_DE_mouse.R -> Performs the DE between Palbociclib treated & 
                                 control samples pseudobulked mouse spots/genes.
                                 
                             INPUT: * data/seurat_rds/all.rds
                                    * data/spot_meta/species_classify_v2/
                                                                   *_species.txt
                            
                             OUTPUT: * data/DE_out/Pseudo_Limma_mouse/*
                           
    X5_clusterProfiler_human.R -> Performs the GSEA on the human DE genes ranked
                                    by t-value.       
                             
                             INPUT: * data/DE_out/Pseudo_Limma_human/
                                               de_results_PseudoLimma_human.xlsx
                             OUTPUT: * data/DE_out/Pseudo_Limma_human/gsea_out/*
                                     * figure_components/DE_figures/
                                                          _human_human_GSEA*.pdf
                                                          
    X6_DE_figures.py -> Creates volcano plots/MDS plots to 
                        visualise the DE results for the human/mouse/mix. 
                        Also includes the PCA analysis.
                        
                        INPUT:  data/DE_out/Pseudo_Limma*/de_results_*
                        
                        OUTPUT: figure_components/DE_figures/*
                                data/DE_out/Table_S2_MixSpotDE.xlsx
                    
    
    
    
    
    
    
    
    
    
    
    
    

    










                   