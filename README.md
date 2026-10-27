
# RoDPAA
RoDPAA (Robust Dissimilarity-based Principal Angle Analysis) is an R package for unsupervised integration of high-dimensional multiview data. It identifies shared sample-level variation across two paired data views through principal-angle analysis of signal subspaces derived from sample dissimilarity matrices.

Unlike methods that require shared variation to occupy an identical subspace across views, RoDPAA allows corresponding latent directions to be partially aligned, recognizing that different data modalities may capture complementary aspects of the same underlying biological processes. By quantifying the strength of cross-view alignment and using permutation-based inference to identify statistically significant shared components, RoDPAA can recover coordinated variation that may be overlooked under an exact common-subspace assumption.

By operating on sample-level dissimilarities rather than the original feature matrices, RoDPAA accommodates user-specified dissimilarity measures, including Euclidean and non-Euclidean measures, and allows each data view to retain a geometry appropriate to its characteristics. This flexibility makes RoDPAA applicable to heterogeneous multiview data, including but not limited to microbiome, metabolomic, transcriptomic, and other high-dimensional biomedical measurements.



## Installation

Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("yanyanxu0428/RoDPAA")
```

## Quick start

RoDPAA needs one dissimilarity matrix per view, computed on the same
participants in the same order. The package includes a small simulated
example: two views measured on 80 participants that share two latent
factors.

``` r
library(RoDPAA)

# 1. one dissimilarity matrix per view
D1 <- dist(rodpaa_example$view1)   # Euclidean dissimilarity, view 1
D2 <- dist(rodpaa_example$view2)   # Euclidean dissimilarity, view 2

# 2. run RoDPAA (the permutation test is random, so set a seed)
set.seed(1)
fit <- rodpaa(D1, D2, input = "dissimilarity")
#> Permutation test selected k = 2
fit
#> RoDPAA fit
#>   participants:             80
#>   signal ranks (r1, r2):    3 , 3
#>   shared dimension (k_hat): 2
#>   leading shared intensities: 0.978, 0.967, 0.081
```

RoDPAA kept three signal directions in each view and found that two of them
are shared. The shared intensities are the cosines of the principal angles
between the two views: values near 1 mean the views vary together across
participants.

``` r
# 3. the permutation test behind k_hat
summary(fit)
#> ... (the fit overview above, then)
#> Permutation test (bonferroni correction, alpha = 0.05):
#>  component singular_value p_value p_adj significant
#>          1      0.9779480   0.000 0.000        TRUE
#>          2      0.9670978   0.000 0.000        TRUE
#>          3      0.0809016   0.163 0.489       FALSE

# 4. see it
plot(fit)
```

<img src="man/figures/README-plot.png" width="70%" />

The first two components lie far above what is expected when participants
are randomly mismatched between the views (grey band), so they are shared.

Two more plots are available: `plot(fit, type = "scores")` shows each
participant's scores on two shared factors (add `group = ` to colour them,
for example by diagnosis), and `plot(fit, type = "pvalue")` shows the
permutation P-values. See `?plot.rodpaa`.

``` r
# 5. the shared factors: one score per participant, for downstream analysis
round(head(fit$shared_factors, 3), 3)
#>        SF1    SF2
#> P01 -0.114 -0.105
#> P02  0.120  0.128
#> P03 -0.123 -0.055
```

## Choosing a dissimilarity

Each view can use its own dissimilarity, for example Bray-Curtis for
microbiome relative abundances and Euclidean distance for log-transformed
metabolite levels. Any dissimilarity can be used, including UniFrac or
Aitchison distance:

``` r
D1 <- vegan::vegdist(microbiome, method = "bray")
D2 <- dist(log(metabolome))
fit <- rodpaa(D1, D2, input = "dissimilarity")
```

**Shortcut: start from the data matrices.** For Euclidean and Bray-Curtis,
RoDPAA can compute the dissimilarities for you with `input = "data"`, choosing
one dissimilarity per view. This runs the same analysis as above:

``` r
fit <- rodpaa(microbiome, log(metabolome), input = "data",
              dissimilarity1 = "bray", dissimilarity2 = "euclidean")
```

RoDPAA uses the data exactly as given. Apply any transformation (for example
CLR or log) and any filtering before calling `rodpaa()`. Bray-Curtis needs
nonnegative data, such as counts or relative abundances.

**Kernels (advanced).** If you already have similarity (kernel) matrices, use
`input = "kernel"`. Kernel matrices must be symmetric. RoDPAA double-centers
them before the analysis, so they need not be centered beforehand. To turn a
dissimilarity matrix into a similarity matrix yourself, use
`dist_to_kernel()`.

## Main options

| Argument | What it does | Default |
|---|---|---|
| `input` | `"dissimilarity"`, `"data"` or `"kernel"` | `"dissimilarity"` |
| `dissimilarity1`, `dissimilarity2` | `"euclidean"` or `"bray"` for each view (data input) | `"euclidean"` |
| `var_thresh` | proportion of variance kept when choosing each signal rank | `0.80` |
| `r1`, `r2` | fix the signal ranks instead of choosing them | chosen from the data |
| `n_perm` | number of permutations | `1000` |
| `alpha` | significance level (Bonferroni-adjusted) | `0.05` |
| `k` | fix the shared dimension and skip the test | tested |
| `k_extract` | extract more directions than are significant, for display | `k_hat` |

See `?rodpaa` for all options and the full output.
