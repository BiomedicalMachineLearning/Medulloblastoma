# Spatial transcriptomics analysis of Medulloblastoma PDX models
For sharing scripts on analyses of 
MedulloBlastoma patient-derived xenograft (MB-PDX) of 
Visium spatial RNA-seq between Palbociclib treated & control samples.

## Installation
Used Python 3.7.10.
Python dependency install:

    pip install -r requirements.txt
    
Need to install stLearn and add to your python path:

     git clone https://github.com/BiomedicalMachineLearning/stLearn.git
     
For R scripts, see session info in docs/r_session.md

## Quick Index

    docs/ -> Contains project documentation.
    
        index.md -> Index which describes each script in scripts/ in terms of 
                    INPUT/OUTPUT, and brief description of function.
                    
    figure_components/ -> Contains the outputted components which are used to 
                            create the paper figures. Is grouped according to 
                            the figure. You can see index.md for which script
                            generates what figure components by checking the 
                            script OUTPUT description. 
                    
    scripts/ -> Contains the scripts described in docs/index.md.
                Scripts are grouped into folders relating to order & analysis
                performed.

        X1_QC_SpeciesClassify/ -> Contains scripts related to spot QC & 
                                           classifying spots by human/mouse/mix.

        X2_DEAnalysis/ -> Contains scripts related to the DE analysis between
                            Palbo treated & control samples pseudobulked by 
                            human/mix/mouse spots. 
                            Also includes PCA & GSEA analyses. 

        X3_GiottoEnrichment/ -> Contains scripts related to the Giotto 
                                    enrichment analysis ('Per-spot enrichment').
                                    
        X4_HumanSpotAnnotation/ -> Contains scripts related to annotating 
                                    human spots
                                    by dominant cell type using reference 
                                    scRNA-seq data using SingleR.
                                    
        X5_MouseSpotAnnotation/ -> Contains scripts related to annotating 
                                    mouse spots
                                    by dominant cell type using reference 
                                    scRNA-seq data using SingleR.
                 

