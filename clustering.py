# script to compute cluster labels for compounds according to their sorption parameters KF and n
# using agglomerative clustering
# Laura Lotteraner, 2023-01-26

import pandas as pd
from sklearn.cluster import AgglomerativeClustering

# import data preparation function
from pmt_sorption_library import data_prep

# define data sources
source_Freundlich = "Tab_1_Freundlich_Parameters.csv"
source_properties = "Tab_2a_Compound_Properties.csv"

original_compounds = pd.read_csv(source_properties, sep = ';', dtype=object)

clusters = pd.DataFrame(columns=['Compound', 'Number of Clusters', 'Cluster Label ACF', 'Cluster Label DeACF', 'Cluster Label OXACF'])

# prepare parameters
sorbents = ['DeACF', 'ACF', 'OXACF']
no_clusters = [2, 3, 4]

# iterate through numbers of clusters
for no_cluster in no_clusters:

    tmp = pd.DataFrame(
        columns=['Compound', 'Number of Clusters', 'Cluster Label ACF', 'Cluster Label DeACF', 'Cluster Label OXACF'])
    tmp['Compound'] = original_compounds['Compound']
    tmp['Number of Clusters'] = no_cluster

    # iterate through all sorbents and compute clustering
    for carbon in sorbents:
        x_sorption = data_prep(source_Freundlich, carbon=carbon)

        agglo = AgglomerativeClustering(n_clusters=no_cluster, linkage='average', affinity = 'euclidean')
        agglo.fit(x_sorption)

        tmp['Cluster Label ' + carbon] = agglo.labels_

    clusters = clusters.append(tmp)

# save cluster labels to csv file
clusters.to_csv('Tab_4_Cluster_Labels_Sorption.csv', sep=';', index = False)