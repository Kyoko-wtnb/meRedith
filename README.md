# meRedith
MEREDITH is an approach to integrate multiple omics data sets (e.g. RNAseq, DNA methylation, miRNA and Copy Number variations) and progect data points into 2D-map by performing [tSNE](https://lvdmaaten.github.io/tsne/).

Please refer [E. Taskesen et al. 2016](http://www.nature.com/articles/srep24949) for detail methods.

## Install meRedith package
The package contains two functions `mederith` and `dbscan_SH`.

1. Install dependecies
..```{r}
install.packages("fpc")
install.packages("cluster")
install.packages("Rtsne")
```
2. Install meRedith
..```(r)
install.packages("devtools")
library(devtools)
install_github("Kyoko-wtnb/meRedith")
```
