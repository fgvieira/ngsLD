# ngsLD

`ngsLD` is a program to estimate pairwise linkage disequilibrium (LD) taking the uncertainty of genotype's assignation into account. It does so by avoiding genotype calling and using genotype likelihoods or posterior probabilities.

### Citation

`ngsLD` is still work under progress...

### Installation

`ngsLD` can be easily installed but has some external dependencies:

* `gcc`: >= 4.9.2 tested on Debian 7.8 (wheezy)
* `zlib`: v1.2.7 tested on Debian 7.8 (wheezy)
* `gsl` : v1.15 tested on Debian 7.8 (wheezy)
* `md5sum`: only needed for `make test`

To install the entire package just download the source code:

    % git clone https://github.com/fgvieira/ngsLD.git

To install these tools just run:

    % cd ngsLD
    % make
    % make test

Executables are built into the main directory. If you wish to clean all binaries and intermediate files:

    % make clean

### Usage

    % ./ngsLD [options] --geno glf/in/file --n_ind INT --n_sites INT --out_prefix output/file

#### Parameters
* `--geno FILE`: input file with genotypes, genotype likelihoods or genotype posterior probabilities.
* `--probs`: is the input genotype probabilities (likelihoods or posteriors)?
* `--log_scale`: is the input in log-scale?
* `--n_ind INT`: sample size (number of individuals).
* `--n_sites INT`: total number of sites.
* `--pos` FILE: input file with site coordinates.
* `--max_kb_dist DOUBLE`: maximum distance between SNPs (in Kb) to calculate LD. If set to 0 (zero) will perform all comparisons. [100]
* `--max_snp_dist INT`: maximum distance between SNPs (in number of SNPs) to calculate LD. If set to 0 (zero) will perform all comparisons. [0]
* `--min_maf DOUBLE`: minimum SNP minor allele frequency. [0.001]
* `--call_geno`: call genotypes before running analyses.
* `--N_thresh DOUBLE`: minimum threshold to consider site; missing data if otherwise (assumes -call_geno).
* `--call_thresh DOUBLE`: minimum threshold to call genotype; left as is if otherwise (assumes -call_geno).
* `--rnd_sample DOUBLE`: proportion of comparisons to randomly sample. [1]
* `--seed INT`: random number generator seed for random sampling (--rnd_sample).
* `--out FILE`: output file name. [stdout]
* `--n_threads INT`: number of threads to use. [1]
* `--version`: prints program version and exits.
* `--verbose INT`: selects verbosity level. [1]

### Input data
As input, `ngsLD` accepts both genotypes, genotype likelihoods (GP) or genotype posterior probabilities (GP). Genotypes must be input as gziped TSV with one row per site and one column per individual (__n_sites\*n_ind__) and genotypes coded as [-1, 0, 1, 2].
As for GL and GP, `ngsLD` accepts both gzipd TSV and binary formats, but with 3 columns per individual (__3\*n_sites\*n_ind__) and, in the case of the binary, the GL/GP coded as doubles

It is advisable that SNPs be called first, since monomorphic sites are not informative and it will greatly speed up computation. If not, these comparisons will show up as `nan` or `inf` in the output.

### Output
`ngsLD` outputs a TSV file with LD results for all pairs of sites for which LD was calculated, where the first two columns are positions of the SNPs, the third column is the distance (in bp) between the SNPs, and the following 4 columns are the various measures of LD calculated (r^2 from pearson correlation between expected genotypes, D from EM algorithm, D' from EM algorithm, and r^2 from EM algorithm).

### Possible analyses
##### LD prunning
For some analyses, linked sites are typically pruned since their presence can bias results. You can use the script `scripts\prune_graph.pl` to prune your dataset and only keep unlinked sites.

##### LD decay
You can also fit an exponential distribution to estimate the rate of LD decay. We provide the script `scripts\Fit_Exp.py` but, for this type of analysis, `--rnd_sample` option should be used since `ngsLD` will be much faster and you don't need all comparisons. The script utilizes the python package `lmfit` (Newville et al., 2016) to run a non-linear least-squares minimisation to fit the model ![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cbig%28%20LD%20%3D%20LD_0%20%5Ccdot%20e%5E%7Bx%20%5Clambda%7D%20%5Cbig%29%20&bc=White&fc=Black&im=png&fs=12&ff=mathdesign&edit=0) to the data; it estimates the coefficients of the exponential decay equation and the Akaike Information Criterion (AIC) for each data set.

    % python Fit_Exp.py --input_type FOLDER --input_name /ngsLD/data/folder/ --data_type r2GLS --plot

* `--input_type STRING`:FILE or FOLDER
* `--input_name STRING`:name of FILE or FOLDER to analyse
* `--data_type STRING`:which measue of LD to analyse. options:r2Pear,r2GLS,D,DPrime. (the two measures of r2 are pearson calculations from expected genotypes and calculations from genotype likelihoods, respectively)
* `--rnd_sample FLOAT`: OPTIONAL argument. Will take a subsample of the available data to use for fitting. Useful for very large files that would otherwise take too long to run. 
* `--plot`: OPTIONAL argument. will use an additional r script to generate a plot of each data set with the fit curve overlaid. 

### Thread pool
The thread pool	implementation was adapted from Mathias Brossard's and is freely available from:
https://github.com/mbrossard/threadpool
