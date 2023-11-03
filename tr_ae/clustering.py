import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import operator
from sklearn.decomposition import PCA
from collections import Counter


class GMM:
    def __init__(self, data, max_k):
        """
        self.X : data to cluster
        self.x,self.y,self.z：first three dimensions of self.X；
        self.max_num_cluster：the maximum number of clusters to consider
        """
        self.X = data
        self.x = self.X[:, 0]
        self.y = self.X[:, 1]
        self.z = self.X[:, 2]
        self.max_num_cluster = max_k

    def cluster(self, path,  cov_type="tied"):
        bic = []
        labels_all = []
        valid = []
        min_bic = []
        n_components = 1
        count = 0
        n_components_range = []
        # Gaussian mixture model is used for clustering, and BIC is used to select the optimal number of clusters
        while n_components <= self.max_num_cluster:
            n_components_range.append(n_components)
            gmm = GaussianMixture(n_components=n_components, covariance_type=cov_type, n_init=100)
            gmm.fit(self.X)
            labels = gmm.predict(self.X)
            labels_u, counts = np.unique(labels, return_counts=True)
            if np.any(counts < 3):
                valid.append(False)
            else:
                valid.append(True)
            labels_all.append([])
            labels_all[n_components - 1].append(labels)
            s = gmm.bic(self.X)
            bic.append(s)
            if n_components == 1 or s < min_bic:
                min_bic = s
                count = 0
            elif s >= min_bic:
                count += 1
            if count >= 20:
                break
            n_components += 1

        n_components_range = np.array(n_components_range)
        n_components_range = n_components_range[valid]
        print(n_components_range)
        bic = np.array(bic)
        bic = bic[valid]
        labels_all = np.array(labels_all)
        labels_all = labels_all[valid,]
        ind, val = min(enumerate(bic), key=operator.itemgetter(1))
        num_cluster = n_components_range[ind]
        labels = np.squeeze(labels_all[ind,])

        print(Counter(labels))
        num_cluster = min(num_cluster, len(Counter(labels)))
        print(num_cluster)

        # visualize clustering results
        plt.figure("3D Scatter", facecolor="lightgray")
        ax3d = plt.gca(projection="3d")
        plt.title('CoT    ', fontsize=20)
        ax3d.set_xlabel('Feature1', fontsize=14)
        ax3d.set_ylabel('Feature2', fontsize=14)
        ax3d.set_zlabel('Feature3', fontsize=14)
        plt.tick_params(labelsize=10)
        ax3d.scatter(self.x, self.y, self.z, s=20, c=[labels], cmap='viridis', marker="o")
        plt.show()

        # The results were displayed by PCA reduced to two dimensions
        pca = PCA(n_components=2)
        pca.fit(self.X)
        X_new = pca.transform(self.X)

        plt.title('CoT                       cluster number:' + str(num_cluster), fontsize=20, loc='left')
        plt.scatter(X_new[:, 0], X_new[:, 1], marker='o', c=[labels])
        plt.xlabel('Feature1', fontsize=20, labelpad=0)
        plt.ylabel('Feature2', fontsize=20, labelpad=0)
        plt.show()

        return labels, num_cluster
