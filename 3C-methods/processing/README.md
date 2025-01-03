### Pipelines
We use 2 different pipelines when processing 3C data depending on the data source:
* If the data is **HiChIP** data, we use [**Juicer**](https://github.com/aidenlab/juicer) to align/filter and then use [**hic2cool**](https://github.com/4dn-dcic/hic2cool) to convert the `.hic` file to the `.cool` format, which is itself just a `.h5` file.
* If the data is **Hi-C** or **Micro-C** data, we use the nextflow pipeline **distiller**, [maintained by the Open2C group](https://github.com/open2c/distiller-nf).
### Normalization
* In the case of **Hi-C**/**Micro-C** data, which is not biased for any particular genomic regions, [ICE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3816492/) (Iterative Correction and Eigenvector Decomposition) is used. \
* For **HiChIP** data, the assumptions used by ICE mean that it cannot be used to normalize the data. Instead, **Square-root Vanilla Coverage** (VC_SQRT) is used, where each element in the matrix is divided by square root of the products of the sums of counts in row $i$, column $j$:
\
$$M_{i,j}=\frac{W_{i,j}}{\sqrt{\sum_{n=0}^N W_{i,n}\times\sum_{n=0}^N W_{n,j}}}$$
