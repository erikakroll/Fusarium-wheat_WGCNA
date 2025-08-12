![Project_Overview](https://github.com/erikakroll/Fusarium-wheat_WGCNA/blob/master/Figure_1_bioinformatics.png)


**Output from dual RNA-seq Weighted Gene Co-expression Network Analysis (WGCNA) studying stage specific expression of Fusarium graminearum infection of wheat**

**Fgram_Module_annotation** 

This includes comma separated value (CSV) files for all genes in each module for both fungal and wheat modules, which are annotated with Module Membership (MM) values, mean FPKM values, 
InterPro annotation and Gene Ontology annotation across stages of infection. 

**Wheat_Module_Annotation** 

This includes comma separated value (CSV) files for all genes in each module for both fungal and wheat modules, which are annotated with Module Membership (MM) values, mean FPKM values, 
InterPro annotation, Plant Trait Ontology (TO) and Gene Ontology annotation across stages of infection. 

**Module_Eigengene** 
The function multiSetMEs from the WGCNA package (Langfelder et al., 2007) was used to calculate module eigengene expression. 

**Network_Construction**
Provides codes and input raw count data for generation of WGCNA networks for F. graminearum and wheat using the WGCNA package (Langfelder et al., 2007). This includes the normalisation of counts. 

**VST_normalised_counts**
CSV files containing VST normalised count data used as input for susquent WGCN analysis. 

**Gene_module_assignments**
List of gene module assignments for all _F graminearum_ genes (ID - RRES v5) and all wheat (_Triticum aestivum_) genes (IDs = IWGSC Refseq v2.1 and v1.1 provided) in all modules. 

**YPD_vs_in_planta_sig_DE**
List of Differentially Expressed (DE) genes from contrasts between PH-1 axenic culture in Yeast Peptone Dextrose (YPD) broth and in planta conditions used in network. 

Langfelder P, Horvath S. 2008. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9:559. doi:10.1186/1471-2105-9-559

RNA-seq data from 

Dilks T, Halsey K, De Vos RP, Hammond-Kosack KE, Brown NA. 2019. Non-canonical fungal G-protein coupled receptors promote Fusarium head blight on wheat. PLoS Pathog 15:e1007666. doi:10.1371/journal.ppat.1007666

available on European Nucleotide Archive: PRJEB75530
