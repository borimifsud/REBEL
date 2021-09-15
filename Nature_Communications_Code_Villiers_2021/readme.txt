# Description of the code detailed in:
# U937_PR9_Processed_Data_Descriptions.R 
# and
# Villiers_et_al_XGBoost.R
# 
# U937_PR9_Processed_Data_Descriptions.R contains a description of each of the U937-PR9 
# datasets used in the study, including the processing (genome alignment, filtering, parameters)
# and how to load the processed datasets for interrogation in R.
# The processed and raw datasets have been deposited on GEO with accession GSE173754
# These are included in the Nature_Communications_Code_Villiers_2021 folder
# U937_PR9_Processed_Data_Descriptions.R can be invoked in R simply by using the following:
# source("./code/U937_PR9_Processed_Data_Descriptions.R") 
# assuming all libraries required are installed
# To install libraries, in R install.packages(library name) can be used
#
# Villiers_et_al_XGBoost.R contains the basic code used for the machine learning predictions
# using the XGBoost in R package.
# The code shows how the motifs were counted on each ATAC-seq fragment using the HOMER tool,
# followed by the initial XGBoost predictions, then how SHAP scores
# were generated for each data point
# The code contains descriptions for each section of the code
# This requires the access to the U937_PR9_Processed_Data_Descriptions.R script and particular
# libraries outlined in the script
#
#
# Overall: The 2 scripts included are to enable:
# 1) the interrogation of each of the NGS datasets generated in this study (U937_PR9_Processed_Data_Descriptions.R)
# 2) to simulate/interrogate the XGBoost experiments performed in this study (Villiers_et_al_XGBoost.R)
# 

