---
aliases:
  - K-means
machineLearningProblem: Clustering
---
## Random initialization trap and k-means++

The k-means++ algorithm solves the random initialization trap. The intuition behind this approach is that spreading out the k initial cluster centers is a good thing. The first cluster center is chosen uniformly at random from the data points that are being clustered, after which each subsequent cluster center is chosen from the remaining data points with probability proportional to its squared distance from the point's closest existing cluster center.

The sklearn implementation of the k-means algorithm provides k-means++ by default.

Automated elbow method: Silhouette, [https://jwcn-eurasipjournals.springeropen.com/articles/10.1186/s13638-021-01910-w](https://jwcn-eurasipjournals.springeropen.com/articles/10.1186/s13638-021-01910-w)