# ngsLD

`ngsLD` is a program to estimate pairwise linkage disequilibrium (LD) taking the uncertainty of genotype's assignation into account. It does so by avoiding genotype calling and using genotype likelihoods or posterior probabilities.

### Citation

`ngsLD` has been submitted to [Bioinformatics](http://bioinformatics.oxfordjournals.org), so please cite it if you use it in your work:

    Fox EA, Wright AE, Fumagalli M, and Vieira FG
    ngsLD: evaluating linkage disequilibrium using genotype likelihoods
    Bioinformatics (2018) ??: ????-????


### Installation

`ngsLD` can be easily installed but has some external dependencies:

* `gcc`: >= 4.9.2 tested on Debian 7.8 (wheezy)
* `zlib`: v1.2.7 tested on Debian 7.8 (wheezy)
* `gsl` : v1.15 tested on Debian 7.8 (wheezy)
* `md5sum`: only needed for `make test`

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
* `--pos` FILE: input file with site coordinates (one per line), where the 1st column stands for the chromosome/contig and the 2nd for the position (bp); remaining columns will be ignored but included in output.
* `--max_kb_dist DOUBLE`: maximum distance between SNPs (in Kb) to calculate LD. If set to 0 (zero) will perform all comparisons. [100]
* `--max_snp_dist INT`: maximum distance between SNPs (in number of SNPs) to calculate LD. If set to 0 (zero) will perform all comparisons. [0]
* `--min_maf DOUBLE`: minimum SNP minor allele frequency. [0.001]
* `--call_geno`: call genotypes before running analyses.
* `--N_thresh DOUBLE`: minimum threshold to consider site; missing data if otherwise (assumes -call_geno).
* `--call_thresh DOUBLE`: minimum threshold to call genotype; left as is if otherwise (assumes -call_geno).
* `--rnd_sample DOUBLE`: proportion of comparisons to randomly sample. [1]
* `--seed INT`: random number generator seed for random sampling (--rnd_sample).
* `--extend_out`: print extended output (see below).
* `--out FILE`: output file name. [stdout]
* `--n_threads INT`: number of threads to use. [1]
* `--verbose INT`: selects verbosity level. [1]

### Input data
As input, `ngsLD` accepts both genotypes, genotype likelihoods (GP) or genotype posterior probabilities (GP). Genotypes must be input as gziped TSV with one row per site and one column per individual ![n_sites.n_ind](http://mathurl.com/ycxtfy8u.png) and genotypes coded as [-1, 0, 1, 2].
As for GL and GP, `ngsLD` accepts both gzipd TSV and binary formats, but with 3 columns per individual ![3.n_sites.n_ind](http://mathurl.com/ycvy5fvx.png) and, in the case of binary, the GL/GP coded as doubles.

It is advisable that SNPs be called first, since monomorphic sites are not informative and it will greatly speed up computation. If not, these comparisons will show up as `nan` or `inf` in the output.

### Output
`ngsLD` outputs a TSV file with LD results for all pairs of sites for which LD was calculated, where the first two columns are positions of the SNPs, the third column is the distance (in bp) between the SNPs, and the following 4 columns are the various measures of LD calculated (![r^2](http://mathurl.com/ya2uo8sp.png) from pearson correlation between expected genotypes, ![D](http://mathurl.com/y8cesmet.png) from EM algorithm, ![D'](http://mathurl.com/y8mgegb8.png) from EM algorithm, and ![r^2](http://mathurl.com/ya2uo8sp.png) from EM algorithm). If the option `--extend_out` is used, then an extra 8 columns are printed with number of samples, minor allele frequency (MAF) of both loci, haplotype frequencies for all four haplotypes, and a chi2 (1 d.f.) for the strength of association ([Collins et al., 1999](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC24792/)).

### Possible analyses
##### LD pruning
For some analyses, linked sites are typically pruned since their presence can bias results. You can use the script `scripts/prune_graph.pl` to prune your dataset and get a list of unlinked sites.

    % perl scripts/prune_graph.pl --in_file testLD_2.ld --max_kb_dist 5 --min_weight 0.5 --out testLD_unlinked.id

* `--in_file FILE`: File with input network [STDIN]
* `--max_kb_dist INT`: Maximum distance between nodes (input file 3rd column) to assume they are connected
* `--min_weight FLOAT`: Minimum weight (in `--weight_field`) of an edge to assume nodes are connected
* `--out FILE`: Path to output file [STDOUT]

For more advanced options, please check script help (`perl scripts/prune_graph.pl --help`).

##### LD decay
If you are interested on the rate of LD decay, you can fit a distribution to your data using the script `scripts/fit_LDdecay.R` to fit LD decay models for ![r^2](http://mathurl.com/ya2uo8sp.png) ([Hill and Weir, 1988](https://www.ncbi.nlm.nih.gov/pubmed/3376052) and [Remington et al., 2001](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC58755/)) and ![D'](http://mathurl.com/y8mgegb8.png) ([Abecassis et al., 2001](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1234912/)).

    % Rscript --vanilla --slave scripts/fit_LDdecay.R --ld_files ld_files.list --out plot.pdf

* `--ld_files FILE`: file with list of LD files to fit and plot (if ommited, can be read from STDIN)
* `--out`: Name of output plot

For more advanced options, please check script help (`Rscript --vanilla --slave scripts/fit_LDdecay.R --help`).

##### LD blocks
To plot LD blocks, we also provide a small script as an example for how it can be easily done in `R` using the `LDheatmap` package (by default, ![r^2](http://mathurl.com/ya2uo8sp.png) is plotted).

    cat testLD_2.ld | bash ../scripts/LD_blocks.sh chrSIM_21 2000 5000

### Hints
* `ngsLD` performance seems to drop considerable under extremely low coverages (<1x); consider these cases only if you have large sample sizes (>100 individuals).
* For some analyses (e.g. LD decay) consider sampling your data (`--rnd_sample`), since `ngsLD` will be much faster and you probably don't need all comparisons.
* For the LD decay, as a rule-of-thumb, consider using at least 10'000 SNPs; check the confidence interval and, if too wide, increase number of SNPs.

### Thread pool
The thread pool	implementation was adapted from Mathias Brossard's and is freely available from:
https://github.com/mbrossard/threadpool
