# Index
Index which describes each script in scripts/ in terms of INPUT/OUTPUT, and 
brief description of function.

X1_QC_SpeciesClassify/

    X1_QC.py -> Performs the QC to filter spots & looks at the distribution of 
                                                               the counts/genes.
             INPUT:  * data/Visium8_{sample_name}_*/*
             OUTPUT: * data/filter_ids/{sample_name}_filtered.txt -> 
                              Contains the barcodes of the spots that passed QC.       
                     * figure_components/QC_figure/*
                     
    X2_generate_seurat_rds.R -> Generates seurat rds with sctransform normalised
                                data which is later used for species scoring &
                                giotto enrichment analysis.
                                NOTE: memory intensive, had to run on a hpc.
                                
                    INPUT: * data/Visium8_{sample_name}_*/*
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
                                    
                     INPUT: * data/Visium8_{sample_name}_*/*
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
                    
X3_GiottoEnrichment/

    X1_giotto_spot_enrich.R -> Using the Giotto package for PAGE per-spot 
                                                            enrichment analysis:
                    https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2#Sec2
                    http://spatialgiotto.rc.fas.harvard.edu/giotto.install.native.html#s_macos
                    
                    To perform enrichment for DE E2F targets & different 
                    MB-SHH transcriptional programs, which were downloaded from
                    Hovestadt, et al. into data/third_party_data/.
                    
                    INPUT: * data/seurat_rds/all.rds
                           * data/Pseudo_Limma_human/gsea_out/
                                                         gsea_results_human.xlsx
                           * data/third_party_data/*_program.txt
                           
                    OUTPUT: * data/giotto_out/*
      
    X2_giotto_figures.py -> Creates the spatial plots of the per-spot enrichment
                          in a consistent way across the samples for comparison.
                          
                    INPUT: * data/Visium8_{sample_name}_*/*
                           * data/filter_ids/*
                           * data/giotto_out/SHH_gsea_v2_scores.txt
                           
                    OUTPUT: * figure_components/giotto_figures/*
                    
    X3_SHH_violins.py -> Creates violin plots showing the shift of the 
                              different cell types in response to the treatment. 
                              
                         INPUT: * data/giotto_out/SHH_gsea_v2_scores.txt
                              
                         OUTPUT: * figure_components/giotto_figures/*
                            
X4_HumanSpotAnnotation/

    X1_SME_normalise.py -> Loads in the data, performs SME normalisation, &
                            saves output.
                            
                            INPUT: * data/Visium8_{sample_name}_*/*
                            OUTPUT: * data/scanpy_h5ads/*_all_species_SME.h5ad

    X2_FetalBrain3_SingleR.R -> Loads in the fetal brain data, & uses this as a
                                reference to label cells by dominant cell type.
                    First download the fetal human brain data to match INPUT:
                    Ref. paper: https://www.nature.com/articles/s41586-020-2157-4
                    Download link: https://db.cngb.org/HCL/gallery.html?tissue=Fetal-Brain3

                    INPUT: * data/scanpy_h5ads/*_all_species_SME.h5ad
                           * data/third_party_data/HCL2020/Fetal-Brain3_dge.txt
                           * data/third_party_data/HCL2020/Fetal-Brain3_Anno.txt
                           
                    OUTPUT: * data/spot_meta/*FetalBrain3singleR_scores.txt
                            * figure_components/HumanAnnot_figures/

    X3_FetalBrain3_umap.py -> Loads in the fetal brain data, & UMAP from seurat 
                    to display the data in a figure.
                  Ref. paper: https://www.nature.com/articles/s41586-020-2157-4
                  Download link: https://db.cngb.org/HCL/gallery.html?tissue=Fetal-Brain3
                  
                  INPUT: * data/third_party_data/HCL2020/Fetal-Brain3_dge.txt
                         * data/third_party_data/HCL2020/Fetal-Brain3_Anno.txt
                         
                  OUTPUT: * data/third_party_data/HCL2020/FetalBrain3.h5ad
                          * figure_components/HumanAnnot_figures/*

    X4_singleR_humanLabel_panels.py -> Label the data according to the human 
                        labels from the Fetal Human 3 reference & output figure.
                      
                      INPUT: * data/scanpy_h5ads/*_all_species_SME.h5ad
                             * data/spot_meta/species_classify_v2/
                                           *FetalBrain3_singleR_scores_human.txt
                                           
                      OUTPUT: * figure_components/HumanAnnot_figures/
                                                  *FetalBrain3Labels_spatial.pdf
    
    
X5_MouseSpotAnnotation/

    X1_format_Vladoiu.py -> Compiles the dev. mouse cerebellum from
                            Vladiou, et al. to use a reference to annotate mouse
                            spots later. Very memory intensive; ran on HPC. 
                         
                            Need to download the supplementary count data from 
                            GEO using the http link and untar so matches the 
                            INPUT description:
                            Ref. Paper: https://www.nature.com/articles/s41586-019-1158-7
                            Data link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118068
                            
                            INPUT: * data/third_party_data/
                                         Vladiou2019_Nature_scRNA/GSE118068_RAW/
                            
                            OUTPUT: * data/third_party_data/
                                            Vladiou2019_Nature_scRNA/
                                                     Vladoiu2019_lateRef.h5ad -> 
                                                             Just P7 cell types.
                                                             
    X2_Vladoiu_SingleR.R -> Performs the automated annotation of the mouse spots
                            using the Vladoiu mouse cerbellum data as a reference.
                            
                            INPUT: * data/scanpy_h5ads/*_all_species_SME.h5ad
                                   * data/third_party_data/
                                              Vladiou2019_Nature_scRNA/
                                                        Vladoiu2019_lateRef.h5ad
                                                        
                            OUTPUT: * data/spot_meta/species_classify_v2/
                                               *Vladoiu_singleR_scores_mouse.txt
                                    * figure_components/MouseAnnot_figures/   
                                    
    X3_singleR_mouseLabel_panels.py -> Label the data according to the 
                                  mouse labels from the Vladiou & output figure.
                          
                          INPUT: * data/scanpy_h5ad/*_all_species_SME.h5ad
                                 * data/spot_meta/species_classify_v2/
                                               *Vladoiu_singleR_scores_mouse.txt
                          OUTPUT: * figure_components/MouseAnnot_figures/
                                                         *VladLabels_spatial.pdf
                                                         
    X4_border_enrichment.py -> Perform FET to test for enrichment of
                        dominant spot cell types for the border spots to provide
                        statistical evidence for astrocytes at the border.

                        INPUT: * data/Visium8_{sample_name}_*/*
                               * data/scanpy_h5ad/*_all_species_SME.h5ad
                               * data/spot_meta/species_classify_v2/
                                               *Vladoiu_singleR_scores_mouse.txt

                        OUTPUT: * data/spot_meta/astro_enrich_stats.xlsx  
                                * figure_components/MouseAnnot_figures/
                                                             *border_spatial.pdf
                        
    X5_gfap_panels.py -> Generating gene plot panels to show GFAP 
                                 expression as prognostic for astrocytes. 
                                 
                       INPUT: data/scanpy_h5ad/*_all_species_SME.h5ad
                              data/spot_meta/species_classify_v2/
                                               *Vladoiu_singleR_scores_mouse.txt
                       OUTPUT: figure_components/MouseAnnot_figures/
                                                               *Gfap_spatial.pdf                                           
    
    


                   