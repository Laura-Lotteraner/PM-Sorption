import pandas as pd
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import numpy as np

from sklearn.cluster import AgglomerativeClustering


def data_prep(source, properties = ['log D', 'mol weight', 'pos charge', 'neg charge', 'water sol', 'arom rings'], carbon = None, sorption = None):
    dataset = pd.read_csv(source, sep = ';', dtype=object)
    dataset = dataset.set_index('Compound')

    if carbon:
        x_cluster = dataset.drop(['7-(Ethylamino)-4-methylcoumarin', 'Carbendazim', 'Pipamperone', 'Tramadol'])
        if sorption == 'KF':
            x_cluster = x_cluster.loc[:, ['KF ' + carbon]].values
            x_cluster = MinMaxScaler().fit_transform(x_cluster)
        else:
            x_cluster = x_cluster.loc[:, ['KF ' + carbon, 'n ' + carbon]].values
            x_cluster = MinMaxScaler().fit_transform(x_cluster)

    else:
        x_cluster = dataset[properties].values
        x_cluster = MinMaxScaler().fit_transform(x_cluster)

    return x_cluster


def distances_clustering(x_cluster, linkage_method, affinity):
    agglo = AgglomerativeClustering(linkage=linkage_method, affinity = affinity, compute_distances=True)
    agglo.fit(x_cluster)

    children = agglo.children_
    children_list = children.tolist()
    distances = agglo.distances_.tolist()

    no_samples = len(x_cluster)

    distance_vector = []

    for i in range(no_samples):
        for j in range(i):

            index_1 = i
            index_2 = j

            criterium = True

            while criterium:
                if ([index_1, index_2] in children_list):
                    distance_vector.append(distances[children_list.index([index_1, index_2])])
                    criterium = False
                elif (([index_2, index_1] in children_list)):
                    distance_vector.append(distances[children_list.index([index_2, index_1])])
                    criterium = False
                else:
                    while ([index_1, index_2] not in children_list) and ([index_2, index_1] not in children_list):
                        for k in range(len(children.tolist())):
                            if index_1 in children[k, :]:
                                tmp1 = no_samples + k
                            if index_2 in children[k, :]:
                                tmp2 = no_samples + k
                        if tmp1 < tmp2:
                            index_1 = tmp1
                        else:
                            index_2 = tmp2

    distance_vector = np.array(distance_vector)

    return distance_vector


def distances_direct(x_cluster):
    distance_vector = []

    for i in range(len(x_cluster)):
        for j in range(i):
            distance_vector.append(np.linalg.norm(x_cluster[i,:] - x_cluster[j,:]))

    distance_vector = np.array(distance_vector)

    return distance_vector
