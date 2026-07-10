import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


def plot(input, labelsels, distance_matrixes):
    rank_names = list(input.keys())
    n_cols = max(1, len(rank_names))

    fig, axes = plt.subplots(1, n_cols, figsize=(4 * n_cols, 3.5), squeeze=False)

    for idx, rank in enumerate(rank_names):
        ax = axes[0, idx]
        distance_matrix = np.asarray(distance_matrixes[rank], dtype=float)
        labels = list(labelsels[rank])

        if len(labels) < 2:
            ax.text(0.5, 0.5, "Not enough samples", ha="center", va="center")
            ax.set_axis_off()
            ax.set_title(f"Rank {rank}")
            continue

        condensed = squareform(distance_matrix, checks=False)
        linkage_matrix = linkage(condensed, method="average")

        dendrogram(
            linkage_matrix,
            ax=ax,
            labels=labels,
            orientation="top",
            leaf_rotation=90,
            leaf_font_size=8,
        )
        ax.set_title(f"Rank {rank}")
        ax.set_xlabel("Samples")
        ax.set_ylabel("Distance")

    fig.suptitle("Dendrogram clustering in different ranks")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
