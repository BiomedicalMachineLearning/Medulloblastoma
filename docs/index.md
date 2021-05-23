# Index
Contains information with regards to the structure of the directory and where
the different files were obtained from. 

info/

    MDB_results_summary_11_03_2021.pptx -> Contains the summary of the results
                                        via powerpoint up to the indicated date.

data/

    Visium*/ -> Obtained from Ryan on the 1st of March, 2021. 
                Contains the Visium data from the PDX models. 
                
    seurat_rds/ -> Contains the seurat objects for each of the spatial data.
    
    ryan_ids/ -> Contains text files for each of the spatial data with the ids
                 after he performed filtering. 
                 Sent at 4:01pm on the 3rd of March via slack.~~~~
                 
    third_party_data/ -> Currently contains the supplementary excel downloaded 
                        from the indicated papers which contain prognostic gene
                        sets. Text files with relevent names are copy-past of 
                        those interesting gene sets so they can be easily read.
         
        Barnes2018_NatureCellBiology_MesynchymalGlioblastoma.xlsx ->
                    Sent by Laura, comes from the paper indicated below.     
                    https://www.nature.com/articles/s41556-018-0183-3#Abs1
                    
        Vladoiu2019_Nature_scRNA/ -> Contains mouse cerebellum developmental
                    scRNA-seq from the indicated paper. Count matrices available
                    at GSE118068.
                    
                    cluster_annotations.csv -> obtained directly from first 
                                            author via email with subject line:
                                            "Childhood cerebellar tumours mirror 
                                       conserved fetal transcriptional programs".
                                       Refers to the plot in extended figure 2.
                                 
    
                                       
    DE_out/ -> Contains the output from running the snakemake pipeline for DE
                (see docs/log.md for reference).
                 
    spot_meta/
            data/spot_meta/cluster_human_mix_vladiouSims.txt ->
                                        Contains the similarity of the pseudo-
                                        bulk of each cluster with the pseudo-
                                        bulk of each of the mouse cerebellum
                                        cell types. (output stlearn_clustering.py).
                                        
            data/spot_meta/clusters_human_mix.txt -> Contains the 
                                        cluster labels and the tissue type 
                                        labels assigned based on the cosine 
                                        similarities of the tissues with 
                                        particular cell types. 
                                        (output stlearn_clustering.py)
                                        
            data/spot_meta/colors.pkl -> Pickled dictionary containing the colors
                                        used when plotting the different stlearn
                                        clusters. (output stlearn_clustering.py)
   
    supps/
        cluster_to_tissues.xlsx -> Contains the cluster to the cell type 
                                    similarities by measuring cosine distance
                                    from vladoiu medulla blastoma cell types. 
                                    (output stlearn_clustering.py)
        cluster_to_tissues_formatted.xlsx -> Same as above but with some excel
                                            formatting. 
scripts/

Analysing data in-&-of itself:
    
    helpers.R -> Contains functions like: loading/saving data, pulling out
                        human/mouse homologues, plotting/ 
    
    seurat_spatial_explore.R -> Went through Seurat v4 spatial vignette in order
                                to understand how to load the data.
                                
    generate_seurat_rds.R -> Goes through each of the spatial datasets, loads
                                them into Seurat, performs SCTransform 
                                normalisation, and then saves to .rds files.
                                Also has a helper function script.
                                OUTPUT: data/seurat_rds/*
                                        all.rds -> main integrated seurat object
                                                    generate on tinaroo. 
                                
    generate_spatial_plots.R -> Generates visualisation of the spatial data,
                                based on the separated data.
                                OUTPUT: figure_components/*_plots/
                                
    merged_gen_spatial_plots.R -> Same as the above, except uses the merged.rds,
                                  which is curently up-to-date with all of the
                                  annotations indicated below.
                                  OUTPUT: figure_components/*_plots/
                                
    human_mouse_classify.R -> Classifies the spots into human, mouse, or mixed
                                based on summing the expression of genes from
                                the different species, plotting these scores,
                                and setting cutoffs. Added to meta data of the
                                Seurat .rds files.
                                OUTPUT: seurat_rds/* -> with updated meta data.
                                        spot_meta/human_mouse_annots.txt -> 
                                                    Contains the human/mouse
                                                    spot annotations.
                                                    
    check_species.R -> Double checks the classification for some spots which 
                        appear outside of the tumour area.  
                                
    species_DE_treatment.R -> Calls DE between treatment and control, 
                                controlling for between replicate variation. 
                                NOTE this script was run on the bioinfr server,
                                not locally!
                                OUTPUT: data/supps/de_results_MDB.xlsx
                                
    broad_GO_analysis.py -> Using the de results from comparing the treatment
                            with the control, now going to perform GSEA in 
                            order to understand the changes.
                            
    giotto_spot_enrich.R -> Using the Giotto package:
                    https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2#Sec2
                    http://spatialgiotto.rc.fas.harvard.edu/giotto.install.native.html#s_macos
                    
                    To perform enrichment for different SHH transcriptional
                    programs, which were downloaded from the nature papers
                    Quan gave me (see QuanPapers/MedullaBlastoma/TuanPapers/).
                    NOTE: The giotto enrichment plots are generated in:
                            merged_gene_spatial_plots.R
                    INPUT: data/third_party_data/
                    OUTPUT: figure_components/giotto_enrich_plots/
                            merged.rds -> Updated to include the enrichment info
                            
    medullablastoma_TI.R -> Goal is to experiment with different ways of 
                    performing trajectory inference in the data, namely 
                    as a way of connecting/showing gradients between the 
                    SHHA/B/C subtypes which I identified using the enrichment
                    analysis above. NOT FINISHED YET

    conditionSpatialDE/ -> Contains scripts for running a snakemake pipeline for
                            DE analysis of spatial data. 
                                
        GSEA_analyse.py -> Performing the GSEA from the results using Limma
                                with pseudo-bulking.  
                            OUTPUT: data/DE_out/TMM_Limma_Voom/gsea_out/*  
                            
        Check_False_Positive.R -> Gets the overlap between DE genes & stably 
                            expressed/HK genes in order to assess FPs.
                            INPUT: data/DE_out/Pseudo_TMM_Limma_Voom/de_results_PseudoLimma.xlsx
                            or     data/supps/de_results_MDB.xlsx 
                            OUTPUT: data/DE_out/Pseudo_TMM_Limma_Voom/
                            or      data/supps/                
    
    Batch_Checks_Stable_Genes.R -> Plotting a PCA of the stable genes, should 
                     show not much difference in the samples if no batch effect.
                     
    stlearn_clustering.py -> Script for obtaining equivalent clusters of the 
                                    data using stlearn with STSME normalisation. 
                           OUTPUT: 
                           data/spot_meta/cluster_human_mix_vladiouSims.txt,
                           data/spot_meta/clusters_human_mix.txt,         
                           data/spot_meta/colors.pkl, cluster_to_tissues.xlsx
              
        stlearn_helpers.py -> Helper functions for using stlearn. 
                                Predominantly helpers stlearn_clustering.py. 
                                
    stlearn_cluster_DE.R -> Script for using the clusters generated from 
                        stlearn_clustering.py in order to get DE markers genes
                        using Seurat. 
                        
    stlearn_umap.py -> Creating a UMAP for visualising the clusters and the 
                                        merged clusters for each of the samples. 
                        
    cluster_TRIAGE.py -> Uses the stlearn clusters, pseudobulks, normalises, &
                            then perform TRIAGE analysis on these to identify 
                            biologically important genes for each cluster to 
                            help justify annotations. 
                         OUTPUT: data/DE_out/stlearn_cluster_DE/
                                                          triage_by_cluster.xlsx
                                                          
    spot_TRIAGE.py -> Using the stlearn clusters, normalises/scales for each 
                        spot & then performs the rank-ordering followed by 
                        FET enrichment for each gene above a particular cutoff
                        within each cluster of spots (comparing frequency about
                        cutoff in list versus below cutoff).
                        (perhaps impliment automatic detection of a good cutoff
                        based on inflection point).a
                        OUTPUT: data/DE_out/stlearn_cluster_DE/
                                                       triage_spot_profiles.xlsx
                                                       
    lr_spot_analysis.py -> Performs the per-spot LR analysis of the visium data
                        for each sample, & subsequently outputs visualisations
                        for the top few LR pairs across the samples !
            
    joint_space.R -> Attempting to use Seurat integration to create a joint
                            space of the stLearn derived clusters in order to 
                            make sure the clusters between samples represent the
                            same thing !  
                     OUTPUT: data/seurat_rds/integrated.rds
                            
    joint_clusters_DE.R -> Taking the integrated data & trying to find marker 
                            genes of those clusters in order to better 
                            characterise them.     
            
    cell_cycle.R -> Scoring cell cycle of the cells in order to help annotate 
                    the clusters.
                    INPUT: data/seurat_rds/integrated.rds
                    
    Cell_Type_Enrichment.R -> Purpose of this script is to perform enrichment 
                    for cell type specific genes per spot in an attempt to label
                    spots in terms of their active pathways. 
                    INPUT: data/seurat_rds/integrated.rds
                    
    Extract_Enrich_Scores.R -> Extracts the Giotto enrichment scores from the 
                        all.rds & saves to a file; this is based on 2nd version
                        of running the Giotto, but will directly save this 
                        scores when I re-run.
                        
            
Scripts for analysing the data with the Vladoiu dataset:       
    
    NOTE: stlearn_clustering.py above actually has the code where measure the 
            similarity of the Valdoiu data with the stlearn clusters 
            (pseudobulked). There is also a script below which perform this 
            spot-wise. 

    format_Vladoiu.py -> Loads in the developing mouse cerebellum scRNA data and
                                 formats to a set of hdf files for fast loading. 
                        *** This was run on the server in the home directory ***
                        NOTE: saving as a pandas hdf not directly compatible for
                            loading in R. Method to do so is presented here:
                            https://stackoverflow.com/questions/44999871/how-can-i-load-a-data-frame-saved-in-pandas-as-an-hdf5-file-in-r
                        
                        INPUT: data/third_part_data/Vladoiu2019_Nature_scRNA/
                                GSE118068_RAW/*
                        OUTPUT: data/third_part_data/Vladoiu2019_Nature_scRNA/
                                counts.hdf & cell_meta.hdf.     
                
    Vladoiu_Create_Seurat.R -> Creates a Seurat object from the Vladiou2019 data,
                        in preparation for sending to the server and performing
                        SCTransform normalisation. 
                        
                        ** Note that this was updated an run on tinaroo **
                        
                        INPUT: data/third_party_data/Vladoiu2019_Nature_scRNA/
                                cluster_annotations.csv & GSE1180688_RAW/*
                        OUTPUT: data/third_party_data/Vladoiu2019_Nature_scRNA/
                                vladoiu.rds
                  
    Vladoiu_Label_Transfer.R -> Performs the label transfer using the 
                                Vladoiu as a reference, and the spatial spots
                                as the target !
                                
                         **** Note that this was deposited on tinaroo to run ***
                         **** The output has been updated to being on tinaroo **
                        INPUT: -> data/third_party_data/Vladoiu2019_Nature_scRNA/
                                vladoiu.rds,
                               -> data/seurat_rds/all.rds
                        OUTPUT: -> data/third_party_data/Vladoiu*/decon_out/* -
                                    has the text files containing the similarity
                                    estimates, can also find on tinaroo.
                                -> data/seurat_rds/all.rds - adds extra assays
                                    to the Seurat object.
                                -> data/seurat_rds/human.rds -> human gene/spot
                                -> data/seurat_rds/mouse.rds -> mouse gene/spot
                                
    Vladoiu_VisTransfer.R -> Visualises the results from the label transfer.
    
    Vladoiu_ExploreLabelTransferCell2LocNNMF.R -> Takes output from cell2LocNNFM.py
                        when running on the cell type estimates from the Seurat
                        label transfer and visualises in a spatial context. 
                        
    Vladoiu_Spot_Sims.py -> Measures the cosine similarity between the human/mix
                            spots & the pseudobulked cell types !
                            OUTPUT: data/spot_meta/vladoiu_spot_sims.txt
                  
Scripts for analysing the data with the Ocasio dataset:
 
    Ocasio_Seurat_Label_Transfer.R -> Performing the label transfer using the 
                        Seurat object provided by Ben 
                        
General scripts:

    cell2LocNNMF.py -> Performs the cell2LocNNMF analysis, was written and 
                        intended to run on tinaroo. Currently only adapted to 
                        Valdoiu data, but will generalise in the future. 

################################################################################
# Update 10th May 2021
Noticed several spots primarily composed of ribosomal genes, going to re-run the
QC analysis & the species determination, along with everything else from the
start in order to alleviate this problem. 

Furthermore, we have additional samples from another mouse model, 
thus warranting a consistent QC approach. 

Starting off this approach with new scripts in the directory:

May_REDO/
    
    QC.py -> Performs the QC looks at the distribution of the counts/genes & 
             also the ribosomal/mt gene captured. 
             OUTPUT: data/spot_meta/{sample_name}_QCd_barcodes.txt -> 
                              Contains the barcodes of the spots that passed QC.
                     data/spot_meta/{sample_name}_QCd_genes.txt -> 
                                 Contains the genes of the spots that passed QC.
                     figure_components/QC_out -> Contains the plots from the QC.
                     
    species_classify.py -> Creates the plots from the species classification. 
            OUTPUT: data/figure_components/species_out/*
            
    DE_figures.R -> Creates the DE figures by re-running the DE but instead of 
                    generating the figures in R, writes out the necessary data
                    which is then read in by DE_figures.py to generate the 
                    figures.
            OUTPUT: figure_components/PseudoLimma/DE_figures/* 
            
    GSEA_figures.R -> Creates the GSEA figure visualisations. 
                    * discarded, ended up going with just the heatmap ! *
            OUTPUT: figure_components/PseudoLimma/DE_figures/*
            
    GSEA_figures.py -> Creates the GSEA figure visualisation.
            OUTPUT: figure_components/PseudoLimma/DE_figures/*
    
                            
    
    
                                
    
    
                                
        
           
           
           
           
           
                   
                   
                   
                             