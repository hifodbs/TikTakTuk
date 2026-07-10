import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

from sklearn.manifold import MDS


def plot(input, labelsels, distance_matrixes):
    rank_names = list(input.keys())
    n_cols = max(1, len(rank_names))

    fig, axes = plt.subplots(1, n_cols, figsize=(4 * n_cols, n_cols), squeeze=False)

    # Unique label mapping setup
    all_labels = []
    for rank in rank_names:
        all_labels.extend(labelsels[rank])
    unique_labels = sorted(list(set(all_labels)))

    for idx, rank in enumerate(rank_names):
        ax = axes[0, idx]
        distance_matrix = distance_matrixes[rank]
        labels = labelsels[rank]

        mapped_labels = [unique_labels.index(lbl) for lbl in labels]

        mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
        X_2d = mds.fit_transform(distance_matrix)

        scatter = ax.scatter(
            X_2d[:, 0], X_2d[:, 1], c=mapped_labels, cmap="tab20", edgecolors="k", s=50
        )
        ax.set_title(f"Rank {rank}")
        ax.set_xlabel("MDS 1")
        ax.set_ylabel("MDS 2")

    fig.suptitle("Clustering in different ranks")

    plt.tight_layout()
    plt.show()