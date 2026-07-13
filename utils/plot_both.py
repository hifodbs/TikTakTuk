import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from sklearn.manifold import MDS


def plot(input, labelsels, distance_matrixes):
    rank_names = list(input.keys())
    n_cols = max(1, len(rank_names))

    fig, axes = plt.subplots(2, n_cols, figsize=(4 * n_cols, 7), squeeze=False)

    all_labels = []
    for rank in rank_names:
        all_labels.extend(labelsels[rank])
    unique_labels = sorted(list(set(all_labels)))

    for idx, rank in enumerate(rank_names):
        ax_scatter = axes[0, idx]
        ax_dendro = axes[1, idx]

        distance_matrix = np.asarray(distance_matrixes[rank], dtype=float)
        labels = list(labelsels[rank])

        n_samples = distance_matrix.shape[0]
        mapped_labels = [unique_labels.index(lbl) for lbl in labels]

        if n_samples == 1:
            # Put the single point right in the center
            X_2d = np.array([[0.0, 0.0]])
        elif n_samples == 2:
            # Put two points separated on the X-axis by their actual Jaccard distance
            dist = distance_matrix[0, 1]
            X_2d = np.array([[-dist/2, 0.0], [dist/2, 0.0]])
        else:
            # Run standard MDS for 3 or more individuals
            mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
            X_2d = mds.fit_transform(distance_matrix)

        unique_points, inverse = np.unique(X_2d, axis=0, return_inverse=True)
        counts = np.bincount(inverse)
        sizes = np.where(counts[inverse] == 1, 50, np.where(counts[inverse] == 2, 1500, 5000))

        ax_scatter.scatter(
            X_2d[:, 0], X_2d[:, 1], c=mapped_labels, cmap="tab20", edgecolors="k", s=sizes
        )
        ax_scatter.set_title(f"Rank {rank}")
        ax_scatter.set_xlabel("MDS 1")
        ax_scatter.set_ylabel("MDS 2")
        ax_scatter.set_aspect("equal", adjustable="box")
        ax_scatter.set_xlim(-1, 1)
        ax_scatter.set_ylim(-1, 1)

        if len(labels) < 2:
            ax_dendro.text(0.5, 0.5, "Not enough samples", ha="center", va="center")
            ax_dendro.set_axis_off()
            ax_dendro.set_title(f"Rank {rank} dendrogram")
            continue

        condensed = squareform(distance_matrix, checks=False)
        linkage_matrix = linkage(condensed, method="average")

        dendrogram(
            linkage_matrix,
            ax=ax_dendro,
            labels=labels,
            orientation="top",
            leaf_rotation=90,
            leaf_font_size=8,
        )
        ax_dendro.set_title(f"Rank {rank} dendrogram")
        ax_dendro.set_xlabel("Samples")
        ax_dendro.set_ylabel("Distance")

    fig.suptitle("Clustering in different ranks")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    plt.show()
