# meRedith
MEREDITH is an approach to integrate multiple omics data sets (e.g. RNAseq, DNA methylation, miRNA and Copy Number variations) and progect data points into 2D-map by performing [tSNE](https://lvdmaaten.github.io/tsne/)<sub>1</sub>.

Please refer [Taskesen et al. 2016](http://www.nature.com/articles/srep24949) for detail methods.

## Install
The package contains two functions `mederith` and `dbscan_SH` with three dependencies. Please install them first if you don't have those libraries. The package also include an example data sets from TCGA <sub>2</sub>. They will be used in the [example](##Example-with-TCGA-data-sets) section.

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

## Example
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
- `samples` : cancer labels of samples

### Preparation of own data
When you perform meredith with your own data, please prepare either data.frame or matrix. All data sets have to have the same number of features. In the example, all data sets have the same number of sample size but the number of features can varies. Most importantly, input data should be normalized propery.

### Run MEREDITH
Since example data sets have samples in column and features in row, data need to be transferred as `Rtsne` uses row as data points. Note that input `data` has to be list object. When you perform MEREDITH for one data set, please provide as list like `data=list(your.data)`. Here we perform `Rtsne` one time as an example, however, we highly recommend to perform at least 100 times to optimize cost function (it might take some time).
```{r}
mer.out <- meredith(data=list(dataGE, dataME, dataCN, dataMIR), transfer = T, nTSNE = 1)
```

### DBSCAN clustering
Perform DBSCAN for 2D-map with optimization of silhouette score.
```{r}
clst.out <- dbscan_SH(mer.out$Y, showplot=TRUE)
```

### Density map
This visualization is not part of the function, however, here are some example to make density map with [ggplot2](http://docs.ggplot2.org/current/#).
```{r}
#Preparation of labels
mer.plot <- data.frame(ID=samples$ID[match(colnames(dataGE), samples$ID)], x=mer.out$Y[,1], y=mer.out$Y[,2],
                       cancerlabel=samples$cancerlabel[match(colnames(dataGE), samples$ID)], cluster=clst.out)
clst.label.plot <- data.frame(cluster=paste("Cluster", unique(mer.plot$cluster[mer.plot$cluster!=0]), sep=""),
                              x=as.numeric(by(mer.plot$x[mer.plot$cluster!=0], mer.plot$cluster[mer.plot$cluster!=0], median)),
                              y=as.numeric(by(mer.plot$y[mer.plot$cluster!=0], mer.plot$cluster[mer.plot$cluster!=0], max)))
cancer.label.plot <- data.frame(cancerlabel=sort(unique(as.character(mer.plot$cancerlabel))),
                                x=as.numeric(by(mer.plot$x, mer.plot$cancerlabel, median)),
                                y=as.numeric(by(mer.plot$y, mer.plot$cancerlabel, median)))

#density plot
library(ggplot2)
# with cluster label
ggplot(mer.plot, aes(x=x,y=y))+stat_density2d(data=mer.plot[mer.plot$cluster!=0,],aes(fill=factor(cluster)), alpha=0.15, geom="polygon", linetype=0, n=100, show.legend = F)+geom_point(aes(color=cancerlabel), alpha=0.8)+geom_text(data=clst.label.plot, aes(x=x,y=y, label=cluster), alpha=0.8, size=5, show.legend=F)+scale_color_manual(values=rainbow(length(unique(mer.plot$cancerlabel))))+theme_bw()+theme(legend.position="none")

# with cancer label
ggplot(mer.plot, aes(x=x,y=y))+stat_density2d(data=mer.plot[mer.plot$cluster!=0,], aes(fill=factor(cluster)), alpha=0.15, geom="polygon", linetype=0, n=100, show.legend = F)+geom_point(aes(color=cancerlabel), alpha=0.8)+geom_text(data=cancer.label.plot, aes(x=x,y=y, label=cancerlabel), alpha=0.8, size=5, show.legend=F)+scale_color_manual(values=rainbow(length(unique(mer.plot$cancerlabel))))+theme_bw()+theme(legend.position="none")

# with both labels
ggplot(mer.plot, aes(x=x,y=y))+stat_density2d(data=mer.plot[mer.plot$cluster!=0,], aes(fill=factor(cluster)), alpha=0.15, geom="polygon", linetype=0, n=100, show.legend = F)+geom_point(aes(color=cancerlabel), alpha=0.5)+geom_text(data=clst.label.plot, aes(x=x,y=y, label=cluster), alpha=0.8, size=5, show.legend=F)+geom_text(data=cancer.label.plot, aes(x=x,y=y, label=cancerlabel), size=4, fontface="bold", show.legend=F)+scale_color_manual(values=rainbow(length(unique(mer.plot$cancerlabel))))+theme_bw()+theme(legend.position="none")

```
Results of full data sets of TGCA can be browsed at [http://pancancer-map.ewi.tudelft.nl/](http://pancancer-map.ewi.tudelft.nl/).

## Citation
Please cite the following airtive when you use MEREDITH.

Taskesen, E., Huisman, S.M.H., Mahfouz, A., Krijthe, J.K., de Ridder, J., van de Stolpe, A., van den Akker, E., Verheagh, W., Reinders, M.J.Y. 2016. Pan-cancer subtyping in a 2D-map shows substructures that are driven by specific combinations of molecular characteristics. *Scientific Reports.* 6:24949. [doi:10.1038/srep24949](https://www.ncbi.nlm.nih.gov/pubmed/27109935).

## Contact
Kyoko Watanabe: k.watanabe@vu.nl
, Erdogan Taskesen: e.taskesen@vu.nl

Dept. Complex Trait Genetics ([CTGlab](http://ctg.cncr.nl/))
, VU University Amsterdam

## References
1. Maaten, L.J.P.v.d. and Hinton, G.E. 2008. Visualizing High-Dimentional Data Using t-SNE. *Journal of Machine Learning Research.* 9, 2579-2605. [link](http://www.cs.toronto.edu/~hinton/absps/tsne.pdf).
2. The Cancer Genome Atlas Research Network, et al. 2013. The Cancer Genome Atlas Pan-Cancer analysis project. *Nat Genet* 45, 1113-1120. [doi:10.1038/ng.2764](https://www.ncbi.nlm.nih.gov/pubmed/24071849).
