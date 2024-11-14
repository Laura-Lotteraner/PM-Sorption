# script to compute the distance correlation between compounds' sorption parameter KF
# and their properties for different combinations of properties
# Laura Lotteraner, 2023-01-26

import pandas as pd
import dcor

# import data preparation function
from pm_sorption_library import data_prep

# define data sources
source_Freundlich = "Tab_1_Freundlich_Parameters.csv"
source_properties = "Tab_2a_Compound_Properties.csv"

# prepare parameters
sorbents = ['DeACF', 'ACF', 'OXACF']
properties = ['log D', 'mol weight', 'pos charge', 'neg charge', 'water sol', 'arom rings', 'V']
property_combinations = []
for i in range(len(properties) + 1):
    for j in range(i):
        property_combinations.append(properties[j:i])

# prepare output
dcor_comparison = pd.DataFrame(columns = ['carbon', 'properties', 'dcor'])

# iterate through sorbent, compute direct distance in KF
for carbon in sorbents:
    x_sorption = data_prep(source_Freundlich, carbon = carbon, sorption = 'KF')

    # iterate through property combinations, compute direct distance in properties
    for props in property_combinations:
        x_properties = data_prep(source_properties, properties=props)

        # compute distance correlation
        dcor_val = dcor.distance_correlation(x_sorption, x_properties)

        dcor_comparison = dcor_comparison.append(
            {'carbon': carbon, 'properties': ', '.join(props), 'dcor': dcor_val}, ignore_index=True)

# save distance correlations to csv file
dcor_comparison.to_csv('Tab_6_Distance_Correlation.csv', sep=';', index = False)
