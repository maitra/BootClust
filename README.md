## BootClust: Bootstrapping for Significance in Clustering

### Authors: Ranjan Maitra and Volodymyr Melnykov

This C package implements the methodology in Maitra, Melnykov and Lahiri (2012) for bootstrapping for assessing significance in the clustering of multidimensional datasets. The procedure compares two models and declares the more complicated model a better candidate if there is significant evidence in its favor.

If you use this package, please cite:

Maitra, R., Melnykov, V. & S. N. Lahiri. (2012) "Bootstrapping for significance of compact clusters in multi-dimensional datasets." *Journal of the American Statistical Association*, 107(497):378-392. DOI: 10.1080/01621459.2011.646935. <https://www.tandfonline.com/doi/abs/10.1080/01621459.2011.646935>

Requirements: `Rmath` (R's mathematical library) and `LAPACK`

To compile, please use at the prompt:

`make`

This will create an executable called `BootClust`

To see usage, type at the prompt

`./BootClust`

#### Usage

`./BootClust  -i <file> -O <dir> -n <int> -p <int> -# <int> -N -k -K -s -q -m -# -l -S`

````
OPTIONS
        -i <file> datafile (with n (row) observations of p (column) attributes
        -O <dir>  output directory
        -n <int>  number of observations (rows)
        -p <int>  number of dimensions (columns)
        -N <int>  number of bootstrap replications (1000 by default)
        -k <int> (minimum) number of clusters under the null hypothesis
        -K <int> (maximum) number of clusters under the alternative hypothesis
        -s        if series of tests should be performed (none by default)
        -q <real> FDR level (for the sequence of tests)
        -m <0|1>  if kmeans or hierarchical clustering is the clustering algorithm
        -# <int>  number of initializations for k-means (if kmeans chosen, default 0 which equates to n *p* k)
        -l <int>  linkage for hierarchical clustering (if chosen) with:
                           1 = Ward (default)
                           2 = single
                           3 = complete
                           4 = average
                           5 = McQuitty
                           6 = median
                           7 = centroid
         -S         should observations in each coordinate be scaled by their ranges (maximum - minimum values) (default: no)
````

#### Output  

Classification vectors and set of produced *p*-values (corresponding files are placed in folder `OUTPUT`)

For details on methodology and examples on performance, please see:

Maitra, R., Melnykov, V. & S. N. Lahiri. (2012) "Bootstrapping for significance of compact clusters in multi-dimensional datasets." *Journal of the American Statistical Association*, 107(497):378-392. DOI: 10.1080/01621459.2011.646935. <https://www.tandfonline.com/doi/abs/10.1080/01621459.2011.646935>
