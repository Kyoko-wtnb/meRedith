# meRedith
MEREDITH is an approach to integrate multiple omics data sets (e.g. RNAseq, DNA methylation, miRNA and Copy Number variations) and progect data points into 2D-map by performing [tSNE](https://lvdmaaten.github.io/tsne/).

Please refer [Taskesen et al. 2016](http://www.nature.com/articles/srep24949) for detail methods.

## Contents
- [Installing meRedith](##install-meredith-package)
- Example with TCGA data sets
- Citation
- References

## Install meRedith package
The package contains two functions `mederith` and `dbscan_SH` with three dependencies. Please install them first if you don't have those libraries. The package also include an example data sets from TCGA. They will be used in the [example](##Example-with-TCGA-data-sets) section.

### Installing dependecies
```{r}
install.packages("fpc")
install.packages("cluster")
install.packages("Rtsne")
library(fpc)
library(cluster)
library(Rtsne)
```
Please make sure that the latest version of Rtsne is installed.

### Installing meRedith
```{r}
install.packages("devtools")
library(devtools)
install_github("Kyoko-wtnb/meRedith")
library(meRedith)
```

## Example with TCGA data sets
### Example data sets
The package contains 5 data sets; 4 data frames for omic data and one for features of samples. The data sets are based on tcga stage 3 used in [Taskesen et al. 2016](http://www.nature.com/articles/srep24949), but only first 500 features are extracted. The number of sample is the same (4434 samples).

To check available data,
```{r}
data(package="meRedith")
```
- `dataGE` : normalized RPKM of RNAseq
- `dataME` : normalized DNA methylation
- `dataMIR` : normalized expression level or microRNA
- `dataCN` : normalized copy number variations
- `sample_labels` : cancer labels of samples

### Preparation of own data
When you perform meredith with your own data, please prepare either data.frame or matrix. All data sets have to have the same number of features. In the example, all data sets have the same number of sample size but the number of features can varies. Most importantly, input data should be normalized propery.

### Run MEREDITH
Since example data sets have samples in column and features in row, data need to be transferred as `Rtsne` uses row as data points. Note that input `data` has to be list object. When you perform MEREDITH for one data set, please provide as list like `data=list(your.data)`. Here we perform `Rtsne` one time as an example, however, we highly recommend to perform at least 100 times to optimize cost function (it might take some time).
```{r}
mer.out <- meredith(data=list(dataGE, dataME, dataCN, dataMIR), transfer = T, nTSNE = 1)
```

### DBSCAN clustering

### Density map


## Citation
Taskesen, E., Huisman, S.M.H., Mahfouz, A., Krijthe, J.K., de Ridder, J., van de Stolpe, A., van den Akker, E., Verheagh, W., Reinders, M.J.Y. 2016. Pan-cancer subtyping in a 2D-map shows substructures that are driven by specific combinations of molecular characteristics. Scientific Reports, 6:24949.

## References
Maaten, L.J.P.v.d. and Hinton, G.E. 2008. Visualizing High-Dimentional Data Using t-SNE. Journal of Machine Learning Research. 9, 2579-2605.
