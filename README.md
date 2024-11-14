# PM-Sorption

Details about the data and the interpretation of the results are in the process of being published. 
This repository is structured as follows:

## Data
- _Tab_1_Freundlich_Parameters_ contains sorption data for 22 different chemical compounds and 3 different sorbents
- _Tab_2a_Compound_Properties_ contains property data for 18 selected compounds
- _Tab_2b_Compound_Properties_Correlation_ correlates these properties to each other
- _Tab_3a_CPCC_Compound_Properties_ contains the cophenetic correlation coefficient (CPCC) for the compound properties
- _Tab_3b_CPCC_Sorption_ contains the cophenetic correlation coefficient (CPCC) for the sorption parameters and each sorbent
- _Tab_4_Cluster_Labels_Sorption_ contains cluster labels for the compounds according to their sorption parameters
- _Tab_5_distance_correlation_vs_correlation_ contains the comparison of direct correlation and distance correlation (Tab_6) between compound properties and the sorption parameter KF
- _Tab_6_Distance_Correlation_ contains distance correlation between compound properties and the sorption parameter KF

## Code
- _pm_sorption_library.py_ provides functions that are used in the other scripts
- _cpcc.py_ calculates Tab_3a from Tab_2a and Tab_3b from Tab_1 using agglomerative clustering with different linkage_methods (similarity measures) and metrics
- _distance_correlation.py_ calculates Tab_6 from Tab_Tab_1 and Tab_2a according to results from Tab_3a
- _clustering.py_ clusters calculates Tab_4 from Tab_1 using agglomerative clustering according to results from Tab_3b
