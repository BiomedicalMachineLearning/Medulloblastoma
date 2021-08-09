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
                             
    
                             
                    
    
                    
    
    
    
    
    
    
    
    
    
    
    
    

    










                   