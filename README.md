# ngsLD

`ngsLD` is a program to estimate pairwise linkage disequilibrium (LD) taking the uncertainty of genotype's assignation into account. It does so by avoiding genotype calling and using genotype likelihoods or posterior probabilities.

### Citation

`ngsLD` has been published in [Bioinformatics](https://doi.org/10.1093/bioinformatics/btz200), so please cite it if you use it in your work:

    Fox EA, Wright AE, Fumagalli M, and Vieira FG
    ngsLD: evaluating linkage disequilibrium using genotype likelihoods
    Bioinformatics (2019) 35(19):3855 - 3856

### Installation

`ngsLD` can be easily installed but has some external dependencies:

* Mandatory:
  * `gcc`: >= 4.9.2 tested on Debian 7.8 (wheezy)
  * `zlib`: v1.2.7 tested on Debian 7.8 (wheezy)
  * `gsl` : v1.15 tested on Debian 7.8 (wheezy)
* Optional (only needed for testing or auxilliary scripts):
  * `md5sum`
  * `Perl` packages: `Getopt::Long`, `Graph::Easy`, `Math::BigFloat`, and `IO::Zlib`
  * `R` packages: `optparse`, `tools`, `ggplot2`, `reshape2`, `plyr`, `gtools`, and `LDheatmap`
  * `python3` packages: `graph-tool`, and `pandas`

To install the entire package just download the source code:

    % git clone https://github.com/fgvieira/ngsLD.git

and run:

    % cd ngsLD
    % make

To run the tests (only if installed through [ngsTools](https://github.com/mfumagalli/ngsTools)):

    % make test

Executables are built into the main directory. If you wish to clean all binaries and intermediate files:

    % make clean

### Usage

    % ./ngsLD [options] --geno glf/in/file --n_ind INT --n_sites INT --out output/file

#### Parameters
* `--geno FILE`: input file with genotypes, genotype likelihoods or genotype posterior probabilities.
* `--probs`: is the input genotype probabilities (likelihoods or posteriors)?
* `--log_scale`: is the input in log-scale?
* `--n_ind INT`: sample size (number of individuals).
* `--n_sites INT`: total number of sites.
* `--pos(H)` FILE: input file with site coordinates (one per line), where the 1st column stands for the chromosome/contig and the 2nd for the position (bp); remaining columns will be ignored but included in output; `--posH` assumes there is a header.
* `--max_kb_dist DOUBLE`: maximum distance between SNPs (in Kb) to calculate LD. Set to `0`(zero) to disable filter. [100]
* `--max_snp_dist INT`: maximum distance between SNPs (in number of SNPs) to calculate LD. Set to `0` (zero) to disable filter. [0]
* `--min_maf DOUBLE`: minimum SNP minor allele frequency. [0.001]
* `--ignore_miss_data`: ignore missing genotype data from analyses (e.g. MAF and haplotype frequency estimation).
* `--call_geno`: call genotypes before running analyses.
* `--N_thresh DOUBLE`: minimum threshold to consider position; if highest GL is lower, set as missing data (assumes -call_geno).
* `--call_thresh DOUBLE`: minimum threshold to call genotype; if highest GL is lower, left as is (assumes -call_geno).
* `--rnd_sample DOUBLE`: proportion of comparisons to randomly sample. [1]
* `--seed INT`: random number generator seed for random sampling (--rnd_sample).
* `--extend_out`: print extended output (see below).
* `--out(H) FILE`: output file name; `--outH` prints header. [stdout]
* `--n_threads INT`: number of threads to use. [1]
* `--verbose INT`: selects verbosity level. [1]

### Input data
As input, `ngsLD` accepts both genotypes, genotype likelihoods (GP) or genotype posterior probabilities (GP). Genotypes must be input as gziped TSV with one row per site and one column per individual ![n_sites.n_ind](http://latex.codecogs.com/png.latex?(n_{sites}\cdot&space;n_{ind})) and genotypes coded as [-1, 0, 1, 2].
As for GL and GP, `ngsLD` accepts both gzipd TSV and binary formats, but with 3 columns per individual ![3.n_sites.n_ind](http://latex.codecogs.com/png.latex?(3\cdot&space;n_{sites}\cdot&space;n_{ind})) and, in the case of binary, the GL/GP coded as doubles.

It is advisable that SNPs be called first, since monomorphic sites are not informative and it will greatly speed up computation. If not, these comparisons will show up as `nan` or `inf` in the output.

### Output
`ngsLD` outputs a TSV file with LD results for all pairs of sites for which LD was calculated, where the first two columns are positions of the SNPs, the third column is the distance (in bp) between the SNPs, and the following 4 columns are the various measures of LD calculated (![r^2](http://latex.codecogs.com/png.latex?r^2) from pearson correlation between expected genotypes, ![D](http://latex.codecogs.com/png.latex?D) from EM algorithm, ![D'](http://latex.codecogs.com/png.latex?D') from EM algorithm, and ![r^2](http://latex.codecogs.com/png.latex?r^2) from EM algorithm). If the option `--extend_out` is used, then an extra 8 columns are printed with number of samples, minor allele frequency (MAF) of both loci, haplotype frequencies for all four haplotypes, and a chi2 (1 d.f.) for the strength of association ([Collins et al., 1999](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC24792/)).

If both `--max_kb_dist` and `--max_snp_dist` are set to `0`, `ngsLD` will output all comparisons, even between different chromosomes/scaffolds/contigs.

### Possible analyses
##### LD pruning
For some analyses, linked sites are typically pruned since their presence can bias results. You can use [prune_graph](https://github.com/fgvieira/prune_graph) to prune your dataset and get a list of unlinked sites. Alternatively, you can also use the auxiliary scripts `scripts/prune_graph.pl` or `scripts/prune_ngsLD.py`, but they are much slower (specially the perl script) and no longer supported.

###### `prune_graph`
```
prune_graph --in testLD_2.ld --weight-field column_7 --weight-filter "column_3 <= 50000 && column_7 >= 0.5" --out testLD_unlinked.pos
```
or, if you have an output with header, you can also do:
```
prune_graph --header --in testLD_2.ld --weight-field "r^2" --weight-filter "dist <= 50000 && r^2 >= 0.5" --out testLD_unlinked.pos
```

For more advanced options, please check help (`prune_graph --help`).


#### LD decay
If you are interested on the rate of LD decay, you can fit a distribution to your data using the script `scripts/fit_LDdecay.R` to fit LD decay models for ![r^2](http://latex.codecogs.com/png.latex?r^2) ([Hill and Weir, 1988](https://www.ncbi.nlm.nih.gov/pubmed/3376052) and [Remington et al., 2001](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC58755/)) and ![D'](http://latex.codecogs.com/png.latex?D') ([Abecassis et al., 2001](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1234912/)) over physical (or genetic) distance.

There are two models implemented for ![r^2](http://latex.codecogs.com/png.latex?r^2) and one for ![D'](http://latex.codecogs.com/png.latex?D') decay. Briefly, the first is derived by adjusting the expected ![r^2](http://latex.codecogs.com/png.latex?r^2) under a drift-recombination equilibrium for finite samples sizes and low level of mutation ([Hill and Weir, 1988](https://www.ncbi.nlm.nih.gov/pubmed/3376052)) with a single parameter (rate of decay):

![E_r^2](http://latex.codecogs.com/png.latex?E\left[r^2\right]=\left[\frac{10&plus;C}{(2&plus;C)(11&plus;C)}\right]\cdot\left[1&plus;\frac{(3&plus;C)(12&plus;12C&plus;C^2)}{n(2&plus;C)(11&plus;C)}\right])

The second formulation is an extension of the [Sved, 1971](https://www.ncbi.nlm.nih.gov/pubmed/5170716) model (expected ![r^2](http://latex.codecogs.com/png.latex?r^2) under drift-recombination equilibrium) with three parameters (rate of decay, maximum observed ![r^2](http://latex.codecogs.com/png.latex?r^2) and minimum observed ![r^2](http://latex.codecogs.com/png.latex?r^2)) to account for the range of observed ![r^2](http://latex.codecogs.com/png.latex?r^2) values:

![E_r^2](http://latex.codecogs.com/png.latex?E\left[r^2\right]=\frac{r^2_{high}-r^2_{low}}{1&plus;C}&plus;r^2_{low})

where the rate of decay ![C](http://latex.codecogs.com/png.latex?C=4N_e\rho).
For ![D'](http://latex.codecogs.com/png.latex?D'), we fit the expectation derived by [Abecassis et al., 2001](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1234912/), assuming a recombination rate of ![cm-Mb](http://latex.codecogs.com/png.latex?1cm=1Mb) and fixing ![D'0_1](http://latex.codecogs.com/png.latex?D'_0=1), resulting also on a three parameter model (rate of decay, maximum ![D'](http://latex.codecogs.com/png.latex?D') and minimum ![D'](http://latex.codecogs.com/png.latex?D')):

![E_D'](http://latex.codecogs.com/png.latex?E\left[D'\right]=D'_{low}&plus;(D'_{high}-D'_{low})\cdot&space;D'_0&space;\cdot&space;(1-\theta)^t)

For more information, please refer to the online published supplementary data.

    % Rscript --vanilla --slave scripts/fit_LDdecay.R --ld_files ld_files.list --out plot.pdf

* `--ld_files FILE`: file with list of LD files to fit and plot (if ommited, can be read from STDIN)
* `--out`: Name of output plot
* `--n_ind`: Only relevant when fitting ![r^2](http://latex.codecogs.com/png.latex?r^2) decay to specify the first model

For more advanced options, please check script help (`Rscript --vanilla --slave scripts/fit_LDdecay.R --help`). The shape of the decay curve can give valuable insight into the sample's biology (e.g.):
* `Intersect` (high): small Ne (or bottleneck), natural selection (local LD)*, non-random mating (or inbreeding), recent admixture
* `Decay Rate` (high): large Ne (or population expansion)*, high mutation rate*, high recombination rate*, random mating (or no inbreeding)*, not recently admixture
* `Asymptote` (high): population structure*, small sample size*

NOTE: Please keep in mind that these are just general trends, and that it all depends on the biology/genetics of the sample. Since LD is influenced by many factors, it is usually less straightforward to derive exact predictions from it.

##### LD blocks
To plot LD blocks, we also provide a small script as an example for how it can be easily done in `R` using the `LDheatmap` package (by default, ![r^2](http://latex.codecogs.com/png.latex?r^2) is plotted).

    % cat testLD_2.ld | bash ../scripts/LD_blocks.sh chrSIM_21 2000 5000

### Hints
* `ngsLD` performance seems to drop considerable under extremely low coverages (<1x); consider these cases only if you have large sample sizes (>100 individuals).
* For some analyses (e.g. LD decay) consider sampling your data (`--rnd_sample`), since `ngsLD` will be much faster and you probably don't need all comparisons.
* For the LD decay, as a rule-of-thumb, consider using at least 10'000 SNPs; check the confidence interval and, if too wide, increase number of SNPs.

### Thread pool
The thread pool implementation was adapted from Mathias Brossard's and is freely available from:
https://github.com/mbrossard/threadpool
