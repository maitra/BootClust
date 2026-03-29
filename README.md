# BootClust: Bootstrapping for Significance in Clustering

### Authors: Ranjan Maitra and Volodymyr Melnykov

Compilation: make

Main program: Bootclust

Important parameters:

`
 Sample size: N
 Dimensionality: p
 Number of resamples: RefDistN
 Minimum number of clusters: minK
 Maximum number of clusters: maxK
 Significance level: alpha
 Linkage for hierarchical clustering: hcrit (1 - Ward, 2 - single, 3 - complete, 4 - average, 5 - McQuitty, 6 - median, 7 - centroid)
 Method of clustering: meth (0 - spherical clusters, use kmeans; 1 - ellipsoidal clusters, use hierarchical clustering)

Output:  classification vectors and set of produced p-values (corresponding files are placed in folder OUTPUT)
`
