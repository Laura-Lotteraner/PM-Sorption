# script to compute the cophenetic correlation coefficient (CPCC) for compounds
# both for sorption parameters KF and n, and for compound properties
# Laura Lotteraner, 2023-01-26

import pandas as pd
from scipy.stats import pearsonr

# import data preparation and distance calculations functions
from pmt_sorption_library import data_prep, distances_clustering, distances_direct

# define data sources
source_Freundlich = "Tab_1_Freundlich_Parameters.csv"
source_properties = "Tab_2a_Compound_Properties.csv"

# prepare parameters
sorbents = ['DeACF', 'ACF', 'OXACF']
similarity_measures = ['average', 'complete', 'single']
metrics = ['euclidean', 'l1']
properties = ['log D', 'mol weight', 'pos charge', 'neg charge', 'water sol', 'arom rings', 'V']
property_combinations = []
for i in range(len(properties) + 1):
    for j in range(i):
        property_combinations.append(properties[j: i])

# prepare output
cpcc_sorption = pd.DataFrame(columns = ['carbon', 'similarity_measure', 'metric', 'cpcc'])
cpcc_properties = pd.DataFrame(columns = ['properties', 'similarity_measure', 'metric', 'cpcc'])

# iterate through similarity measures, metrics
for sim_meas in similarity_measures:
    for metric in metrics:
        # iterate through sorbents, compute CPCC for sorption
        for carbon in sorbents:
            x_sorption = data_prep(source_Freundlich, carbon = carbon)

            distances_sorption_clustering = distances_clustering(x_sorption, linkage_method = sim_meas, affinity = metric)
            distances_sorption_direct = distances_direct(x_sorption)

            cpcc, _ = pearsonr(distances_sorption_clustering, distances_sorption_direct)

            cpcc_sorption = cpcc_sorption.append({'carbon': carbon, 'similarity_measure': sim_meas, 'metric': metric, 'cpcc': cpcc}, ignore_index=True)

        # iterate through property combinations, compute CPCC for properties
        for props in property_combinations:
            x_properties = data_prep(source_properties, properties = props)

            distances_properties_clustering = distances_clustering(x_properties, linkage_method=sim_meas, affinity=metric)
            distances_properties_direct = distances_direct(x_properties)

            cpcc, _ = pearsonr(distances_properties_clustering, distances_properties_direct)

            cpcc_properties = cpcc_properties.append(
                {'properties': ', '.join(props), 'similarity_measure': sim_meas, 'metric': metric, 'cpcc': cpcc}, ignore_index=True)

# save CPCC to csv files
cpcc_properties.to_csv('Tab_3a_CPCC_Properties.csv', sep = ';', index = False)
cpcc_sorption.to_csv('Tab_3b_CPCC_Sorption.csv', sep = ';', index = False)






