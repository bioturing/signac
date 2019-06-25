# Signac - A versatile package for Single-cell RNA-Seq analysis

We introduce **Signac**, a versatile R package to facilitate the analysis workflow for single-cell data. It helps to find marker genes faster and more accurate, search for cells with similar expression profiles, integrate multiple datasets in the BioTuring Browser database ([know more about BioTuring Browser](https://bioturing.com/product/bbrowser)), etc. For users with a limited computational resource, we provide the helper functions to exercise all analyses for the large-scale datasets from disk. Because of its speed and flexibility, it can be adapted to any existing R analysis pipeline to help explore single-cell data more efficient.

This package can also be used as the reference for the computational methods used in **BioTuring Browser** software.

### 24/06/2019

Please visit [this blog](https://blog.bioturing.com/) to read about **Venice**, the fast and accurate method for finding marker genes, which is incorporated into Signac. Venice function can be accessed via `Signac::VeniceMarker`

Below is the benchmark when finding 4 types of DE genes in a simulated dataset (using scDD R package). Total 15 methods included Venice are tested. The dataset has 2000 DE genes that are even divided into 4 groups (DE, DP, DB, DM), and 18000 non-DE genes (EP and EE).

  <img src="https://lh3.googleusercontent.com/_1GzN-XnIUNepbLaaJHFQjKMwAssfEOT4loEw-SARGu98qvgEWA1BslCNQb3daM6W6j1_Nn3pg6XqvQAxZOl4ACJuYQuUgPmlBb778_DDrgwIsH1VWtkqyJ-2h2S0vy_6NnWghEO_7TTMRp-YCpTQSTDC9gbQQJa6rGRGhbg8sJiqBXYPgPeyzTVDhe69sBjYB6HxHX6kCKI3HL9Vo9hPRENPtvsa9rGkFx0Ud0OgMxKoqSFArymODACMeZzkVjBX8hYDEcCkRga7JYNXREGCHVRsNSwZ9j-5Nc7yGU_VmxUo-PNTd-VFu7z0upcJyvkx4013gBtmH0BncB0dQfI6fVMks8OCBCQuxhojVVLZK0nV7FwqiPR645gdHk9Os8driY50jFgxJuwCtNabAR33UC_QzAB0Pc2ah7dYeAwNfX7dtfCiS-OsnVthNnJxBxaDL2uBBXRXLXXhUt7S1AC7iHw5cmr-4bW5XHUg_1aQKeAyNy0c_CtvUg7pt6lyGREAAKLtGKFDFfkQgL3bm_5dCVwX518kX4GmHgFq_kGMi5WWuGXE55yxWk1P_Q4TY01CbrpW7mZfaRcXRYtcuW32vIPhVqI20e3PMo55gRoA6O_2rjozK6w2fDI1JO7tjYvU72tJfrcoUprTb3lD4K0evDEXpOPV_Q=w1508-h953-no-tmp.jpg" width="800" style="text-align:center"/> 

## Installation

```R
devtools::install_github("bioturing/signac")
```

## Usage
### Find marker genes with ```Signac::VeniceMarker```
```R
> class(pbmc.mat) 
[1] "dgCMatrix"
attr(,"package")
[1] "Matrix" 
> head(clusters) 
[1] "Memory CD4 T" "B"            "Memory CD4 T" "CD14+ Mono"   "NK"           "Memory CD4 T" 
### Find markers for "NK" cluster
> markers <- Signac::VeniceMarker(pbmc.mat, cluster = (clusters != "NK") + 1)
> head(markers) # Type 1: Up-regulated; -1: Down-regulated
 Gene.Name Similarity Log.P.value P.adjusted.value      Type
1      NKG7  0.1366077   -552.6917    1.277099e-236 0.9870968
2      GNLY  0.1575343   -520.4139    6.656220e-223 1.0000000
3      GZMB  0.1730481   -502.0630    4.138385e-215 1.0000000
4      PRF1  0.2285456   -418.0301    9.703964e-179 1.0000000
5      CST7  0.2894125   -346.8274    6.500973e-148 1.0000000
6      CTSW  0.2838700   -345.5459    1.951351e-147 0.9869281
```
