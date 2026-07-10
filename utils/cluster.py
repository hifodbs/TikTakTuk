from sklearn.cluster import AgglomerativeClustering
import numpy as np
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import cdist




def get_number_clusters_from_silhouette(max_k, distance_matrix, linkage):
    silhouette_scores = []
    cluster_range = range(2, max_k)

    for k in cluster_range:
        hierarchical_cluster = AgglomerativeClustering(
            n_clusters=k, metric="precomputed", linkage=linkage
        )
        labels = hierarchical_cluster.fit_predict(distance_matrix)
        score = silhouette_score(distance_matrix, labels, metric="precomputed")
        silhouette_scores.append(score)

    k_best = cluster_range[np.argmax(silhouette_scores)]
    return k_best


def get_best_k(max_k, distance_matrix, linkage,type="silhouette"):
    #it's possible to add other alogithrms to find out best number of clusters
    if type == "silhouette":
        return get_number_clusters_from_silhouette(max_k, distance_matrix, linkage)
    else:
        return get_number_clusters_from_silhouette(max_k, distance_matrix, linkage)
    

def clusterize_patients(input, distance_metric, linkage, best_k_metric, weights):
    distance_matrixes = {}
    labelsels = {}
    for rank in input.keys():
        X = np.array(list(input[rank].values()))
        if len(X) == 1:
            labelsels[rank] = [1]
            distance_matrixes[rank] = np.zeros((1,1))
            continue
        elif len(X) == 2:
            distance_matrixes[rank] = cdist(X, X, metric=distance_metric,w=weights)
            if distance_matrixes[rank][0][1] < 0.5:
                labelsels[rank] =[1, 1]
            else:
                labelsels[rank] =[1, 2]
            continue
        distance_matrix = cdist(X, X, metric=distance_metric,w=weights)

        k = get_best_k(len(X), distance_matrix, linkage, best_k_metric)
        print("foudn best k",k)
        hierarchical_cluster = AgglomerativeClustering(
            n_clusters=k, metric="precomputed", linkage=linkage
        )

        labelsels[rank] = hierarchical_cluster.fit_predict(distance_matrix)
        distance_matrixes[rank] = distance_matrix
    return distance_matrixes, labelsels
