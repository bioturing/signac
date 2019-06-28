---
layout: page 
title: "Tutorial 1: Use Venice to find marker genes"
subtitle: "Simple usage to find marker genes for all clusters in a Single-cell data"
image: /assets/images/Tutorial_1__Use_Venice_to_find_marker_genes_files/figure-markdown_strict/unnamed-chunk-6-1.png
permalink: /Use-Venice-to-find-marker-genes/
---
In this vignette, we will find marker genes for every cluster in the
**PBMC3k** dataset. For convenience, we use the pre-processed R object
from Seurat.v3 tutorial (<a href="https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html" class="uri">Guided Clustering</a>)
to find markers for the annotated clusters. You can download the object
[here](https://www.dropbox.com/s/63gnlw45jf7cje8/pbmc3k_final.rds?dl=1).
This vignette requires Signac and Seurat.v3 has been installed. If you
have not installed **Signac**, please install the latest version at
<a href="https://github.com/bioturing/signac" class="uri">https://github.com/bioturing/signac</a>

```R
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(require(dplyr))
```

First let take a look at the structure of the pre-processed object

```R
pbmc <- readRDS("~/Downloads/pbmc3k_final.rds")
pbmc

## An object of class Seurat 
## 13714 features across 2638 samples within 1 assay 
## Active assay: RNA (13714 features)
##  2 dimensional reductions calculated: pca, umap

table(Idents(pbmc))

## 
##  Naive CD4 T Memory CD4 T   CD14+ Mono            B        CD8 T 
##          697          483          480          344          271 
## FCGR3A+ Mono           NK           DC           Mk 
##          162          155           32           14
```

Now we find marker genes for all annotated clusters using `VeniceAllMarkers` function in Signac.

```R
system.time(pbmc.markers <- VeniceAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 0.25, verbose = F))

##    user  system elapsed 
##   5.258   0.279   5.540

pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = Log2.fold.change)

## # A tibble: 18 x 8
## # Groups:   cluster [9]
##    Gene.Name Log2.fold.change Log10.adjusted.… Log10.p.value cluster
##    <fct>                <dbl>            <dbl>         <dbl> <fct>  
##  1 CD3D                 0.631           -55.6          -58.0 Naive …
##  2 CCR7                 0.639           -26.9          -28.9 Naive …
##  3 IL32                 0.738           -83.7          -87.9 Memory…
##  4 IL7R                 0.666           -35.0          -37.6 Memory…
##  5 S100A9               1.95           -288.          -292.  CD14+ …
##  6 S100A8               2.00           -281.          -285.  CD14+ …
##  7 CD79A                1.83           -210.          -213.  B      
##  8 CD79B                1.54           -178.          -181.  B      
##  9 CCL5                 1.57           -140.          -144.  CD8 T  
## 10 NKG7                 1.50           -137.          -141.  CD8 T  
## 11 LST1                 1.43            -94.0          -98.2 FCGR3A…
## 12 FCGR3A               1.59            -82.8          -86.2 FCGR3A…
## 13 GZMB                 2.00           -103.          -107.  NK     
## 14 GNLY                 2.06           -103.          -107.  NK     
## 15 HLA-DQA1             1.40            -12.4          -16.5 DC     
## 16 FCER1A               1.60            -11.6          -15.2 DC     
## 17 PPBP                 2.72             -9.11         -13.2 Mk     
## 18 PF4                  2.54             -9.11         -13.2 Mk     
## # … with 3 more variables: Dissimilarity <dbl>, Bin.count <dbl>,
## #   Up.Down.score <dbl>
```

Compare with `Seurat::FindAllMarkers` results

```R
system.time(pbmc.markers.seurat <- FindAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 0.25, verbose = F, min.cells.feature = 0))

##    user  system elapsed 
## 114.057   1.416 115.483

pbmc.markers.seurat %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

## # A tibble: 18 x 7
## # Groups:   cluster [9]
##        p_val avg_logFC pct.1 pct.2 p_val_adj cluster      gene    
##        <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>        <chr>   
##  1 1.61e- 82     0.922 0.436 0.11  2.20e- 78 Naive CD4 T  CCR7    
##  2 7.37e- 30     0.764 0.245 0.084 1.01e- 25 Naive CD4 T  LDLRAP1 
##  3 7.95e- 89     0.892 0.981 0.642 1.09e- 84 Memory CD4 T LTB     
##  4 1.85e- 60     0.859 0.422 0.11  2.54e- 56 Memory CD4 T AQP3    
##  5 0.            3.86  0.996 0.215 0.        CD14+ Mono   S100A9  
##  6 0.            3.80  0.975 0.121 0.        CD14+ Mono   S100A8  
##  7 0.            2.99  0.936 0.041 0.        B            CD79A   
##  8 9.48e-271     2.49  0.622 0.022 1.30e-266 B            TCL1A   
##  9 2.96e-189     2.12  0.985 0.24  4.06e-185 CD8 T        CCL5    
## 10 2.57e-158     2.05  0.587 0.059 3.52e-154 CD8 T        GZMK    
## 11 3.51e-184     2.30  0.975 0.134 4.82e-180 FCGR3A+ Mono FCGR3A  
## 12 2.03e-125     2.14  1     0.315 2.78e-121 FCGR3A+ Mono LST1    
## 13 7.95e-269     3.35  0.961 0.068 1.09e-264 NK           GZMB    
## 14 3.13e-191     3.69  0.961 0.131 4.30e-187 NK           GNLY    
## 15 1.48e-220     2.68  0.812 0.011 2.03e-216 DC           FCER1A  
## 16 1.67e- 21     1.99  1     0.513 2.28e- 17 DC           HLA-DPB1
## 17 7.73e-200     5.02  1     0.01  1.06e-195 Mk           PF4     
## 18 3.68e-110     5.94  1     0.024 5.05e-106 Mk           PPBP

```

Finally, we use `DoHeatmap` function from **Seurat** package to draw two
heatmaps of expression of the marker genes found by two method: Seurat
default and Harmony to see the distinct expression pattern of each cell
type (cluster). We only plot top 20 features (all features if less than
20).

### Heatmap of marker genes found by Venice

```R
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = -Log10.adjusted.p.value)
DoHeatmap(pbmc, features = as.character(top10$Gene.Name)) + NoLegend()
```
![Marker genes found by Venice](/assets/images/Tutorial_1__Use_Venice_to_find_marker_genes_files/figure-markdown_strict/unnamed-chunk-6-1.png)

### Heatmap of marker genes found by Seurat default method
```R
top10 <- pbmc.markers.seurat %>% group_by(cluster) %>% top_n(n = 10, wt = -p_val_adj)
DoHeatmap(pbmc, features = as.character(top10$gene)) + NoLegend()
```
![Marker genes found by Seurat default method](/assets/images/Tutorial_1__Use_Venice_to_find_marker_genes_files/figure-markdown_strict/unnamed-chunk-7-1.png)

To perform the Differential expression (DE) test on two individual
clusters, one can use the function `VeniceFindMarkers`. We designed the
user interface that is similar to Seurat R package (Butler et al.,
Nature Biotechnology 2018). Please specify two cluster names `ident.1`,
`ident.2` as the parameters.

### Find DE genes between CD14+ Mono and FCGR3A+ Mono using Venice algorithm
```R
head(Signac::VeniceFindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", logfc.threshold = log(2)))

##   Gene.Name Log2.fold.change Log10.adjusted.p.value Log10.p.value
## 2    FCGR3A       -1.6631222              -71.17833     -75.01446
## 3    S100A8        1.2752262              -62.41254     -66.07258
## 4    S100A9        0.7930380              -61.23949     -64.77459
## 5    IFITM2       -0.8010430              -56.30591     -59.74410
## 6    LGALS2        1.2904962              -49.87441     -53.23342
## 7      GPX1        0.7329328              -39.49077     -42.78284
##   Dissimilarity Bin.count Up.Down.score
## 2     0.7303167         3            -1
## 3     0.6406856         5             1
## 4     0.6282687         5             1
## 5     0.5724071         4            -1
## 6     0.5219134         3             1
## 7     0.4105139         4             1
```
### Find DE genes between CD14+ Mono and FCGR3A+ Mono using Seurat default algorithm
```R
    head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", logfc.threshold = log(2)))

    ##                p_val avg_logFC pct.1 pct.2    p_val_adj
    ## FCGR3A 1.193617e-101 -2.617707 0.131 0.975 1.636926e-97
    ## LYZ     8.134552e-75  1.812078 1.000 0.988 1.115572e-70
    ## RHOC    4.479768e-68 -1.611576 0.162 0.864 6.143554e-64
    ## S100A8  7.471811e-65  2.610696 0.975 0.500 1.024684e-60
    ## S100A9  1.318422e-64  2.286734 0.996 0.870 1.808084e-60
    ## IFITM2  4.821669e-64 -1.445772 0.677 1.000 6.612437e-60
```
