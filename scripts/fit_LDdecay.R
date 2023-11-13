#!/bin/env Rscript

#FileName: fit_LDdecay.R v1.1.0
#Author: "Filipe G. Vieira (fgvieira _at_ sund [dot] ku [dot] dk)"
#Author: "Emma Fox (e.fox16 _at_ imperial [dot] ac [dot] uk)"

##IMPORTS
library(optparse)
library(lgr)
library(jsonlite)
library(tools)
suppressMessages(library(ggplot2))
library(reshape2)
library(plyr)
options(width = 200) 



##SPECIFY ARGUMENTS
option_list <- list(
  make_option(c('--ld_files'), action='store', type='character', default=NULL, help = 'File with list of LD files to fit and plot (if ommited, can be read from STDIN)'),
  make_option(c('--header'), action='store_true', type='logical', default=FALSE, help='Input file has header'),
  make_option(c('--col'), action='store', type='numeric', default=3, help='Which column is distance between sites? [%default]'),
  make_option(c('--ld'), action='store', type='character', default="r2", help='Which LD stats to plot (r2_ExpG, D, Dp, r2) [%default]'),
  make_option(c('--n_ind'), action='store', type='numeric', default=0, help='Number of individuals (for 1-parameter r^2 fitting correction)?'),
  make_option(c('-r', '--use_recomb_rate'), action='store_true', type='logical', default=FALSE, help='Assume constant recombination rate. [%default]'),
  make_option(c('--recomb_rate'), action='store', type='numeric', default=1, help='Recombination rate (or probability of recombination between adjacent sites in cM/Mb) to calculate genetic distances from physical distances. It is assumed to be constant throughout the whole dataset and, for human datasets, a common rule-of-thumb value is 1cM/Mb (1e-6). [%default]'),
  make_option(c('--min_maf'), action='store', type='numeric', default=0, help='Minimum MAF of either site to include in the fitting analysis. [%default]'),
  make_option(c('--max_kb_dist'), action='store', type='numeric', default=Inf, help='Maximum distance between SNPs (in kb) to include in the fitting analysis. [%default]'),
  make_option(c('--fit_boot'), action='store', type='numeric', default=0, help='Number of bootstrap replicates for fitting CI. [%default]'),
  make_option(c('--fit_bin_size'), action='store', type='numeric', default=250, help='Bin data into fixed-sized windows for fitting. [default %default bps]'),
  make_option(c('--fit_level'), action='store', type='numeric', default=1, help='Fitting level: no fitting (0), BFGS (1), best of BFGS/Nelder-Mead (2), and best of BFGS/Nelder-Mead/L-BFGS-B (3). [%default]'),
  make_option(c('--plot_group'), action='store', type='character', default='File', help='Group variable'),
  make_option(c('--plot_data'), action='store_true', type='logical', default=FALSE, help='Also plot data points?'),
  make_option(c('--plot_bin_size'), action='store', type='numeric', default=0, help='Bin data into fixed-sized windows for plotting. [default %default bps]'),
  make_option(c('--plot_x_lim'), action='store', type='numeric', default=NULL, help='X-axis plot limit (in kb). [%default]'),
  make_option(c('--plot_y_lim'), action='store', type='numeric', default=NULL, help='Y-axis plot limit. [%default]'),
  make_option(c('--plot_axis_scales'), action='store', type='character', default='fixed', help='Plot axis scales: fixed (default), free, free_x or free_y'),
  make_option(c('--plot_size'), action='store', type='character', default='1,2', help='Plot size (height,width). [%default]'),
  make_option(c('--plot_scale'), action='store', type='numeric', default=1.5, help='Plot scale. [%default]'),
  make_option(c('--plot_wrap'), action='store', type='numeric', default=0, help='Plot in WRAP with X columns (default in GRID)'),
  make_option(c('--plot_no_legend'), action='store_true', type='logical', default=FALSE, help='Remove legend from plot'),
  make_option(c('--plot_shapes'), action='store_true', type='logical', default=FALSE, help='Use also shapes (apart from colors)'),
  make_option(c('--plot_line_smooth'), action='store', type='numeric', default=1000, help='LD decay curve smoothness'),
  make_option(c('--bin_quant'), action='store', type='numeric', default=0, help='Quantile to represent the bins (e.g. 0 = mean, 50 = median). [%default]'),
  make_option(c('-f','--plot_wrap_formula'), action='store', type='character', default=NULL, help='Plot formula for WRAP. [%default]'),
  make_option(c('-o','--out'), action='store', type='character', default=NULL, help='Output file'),
  make_option(c('--seed'), action='store', type='numeric', default=NULL, help='Seed for random number generator'),
  make_option(c('-l','--loglevel'), action="store", type="character", default="info", help='Log verbosity level. [%default]')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Set log
lgr::threshold(opt$loglevel)

# Set random seed
if(is.null(opt$seed))
  opt$seed <- as.integer(runif(1,1,100000))
lgr$info("Random seed: %i", opt$seed)
set.seed(opt$seed)

# Parse input LD files
if(is.null(opt$ld_files)) 
  opt$ld_files <- file('stdin')
ld_files <- read.table(opt$ld_files, header=opt$header, stringsAsFactors=FALSE)
colnames(ld_files)[1] <- "File"
# Keep file list ordered
for (id in colnames(ld_files))
  ld_files[,id] <- factor(ld_files[,id], unique(ld_files[,id]), ordered=TRUE)
n_files <- nrow(ld_files)
lgr$debug("Input files", data = toJSON(ld_files))

# Parse LD stats to plot
if(!is.null(opt$ld))
  opt$ld <- unlist(strsplit(opt$ld, ","))
if(!all(opt$ld %in% c("r2_ExpG", "D", "Dp", "r2")))
  stop(lgr$error("Invalid LD measure to plot"))
n_ld = length(opt$ld)

# Check if number of individuals was specified
if(opt$n_ind < 0)
  stop(lgr$error("Number of individuals must be greater than zero"))
if(opt$n_ind > 0 && length(unique(ld_files$File)) > 1)
  lgr$warn("Sample size for correction will be assumed equal for all samples!")
if(!any(opt$ld %in% c("r2","r2_ExpG")) && opt$n_ind)
  stop(lgr$error("Number of individuals is only used for r^2 fitting"))
for(i in opt$ld)
  lgr$info("Fitting %s LD decay assuming a %s parameter decay model", i, ifelse(i == "Dp" || opt$n_ind == 0, "three (rate of decay, max LD and min LD)", "one (rate of decay)"))

# Check max_kb_dist parameter 
if(opt$max_kb_dist < 50)
  lgr$warn("Fitting of LD decay is highly unreliable at short distances (<50kb).")

# Check binning sizes
if(opt$plot_bin_size > 0 && opt$plot_bin_size < opt$fit_bin_size)
   stop(lgr$error("Ploting bin size must be greater than fiting bin size!"))

# Parse plot_size
opt$plot_size <- as.numeric(unlist(strsplit(opt$plot_size, ",")))
if(length(opt$plot_size) < 2)
  opt$plot_size <- c(opt$plot_size, opt$plot_size)

# Plot formula
if(!is.null(opt$plot_wrap_formula))
  opt$plot_wrap_formula <- as.formula(opt$plot_wrap_formula)

# Set output file name (if not defined)
if(is.null(opt$out)){
  if(is.null(opt$ld_files))
    stop(lgr$error('Output file name required, when reading LD files from STDIN'))
  opt$out <- paste(basename(file_path_sans_ext(opt$ld_files)),".pdf", sep = "")
}


### Load LD data
ld_data <- data.frame()
for (i in 1:n_files) {
  ld_file <- as.character(ld_files$File[i])
  # Read point data
  lgr$info("Reading file %s", ld_file)
  tmp_data <- read.table(gzfile(ld_file), sep="\t", header=TRUE, quote="\"", dec=".")
  lgr$info("Number of lines: %i", nrow(tmp_data))
  lgr$debug("Data", data = toJSON(head(tmp_data)))
  # Check if file is valid
  if(ncol(tmp_data) < 5)
    stop(lgr$error('Invalid LD file format.'))
  # Filter by minimum MAF
  if ("maf1" %in% colnames(tmp_data) & "maf2" %in% colnames(tmp_data)) {
    tmp_data <- tmp_data[which(tmp_data$maf1 >= opt$min_maf & tmp_data$maf2 >= opt$min_maf),]
  }
  # Filter by maximum distance
  tmp_data <- tmp_data[which(tmp_data$dist < opt$max_kb_dist*1000),]
  lgr$info("Number of lines after filtering: %i", nrow(tmp_data))
  # Extract relevant columns
  tmp_data <- tmp_data[, which(names(tmp_data) %in% c("dist",opt$ld))]
  # Convert all 'Inf' to NA
  tmp_data[mapply(is.infinite, tmp_data)] <- NA
  # Calculate genetic distances, according to Haldane's formula (assumes constant rate across all dataset)
  if(opt$use_recomb_rate && !is.null(opt$recomb_rate))
    tmp_data$dist <- (1 - (1 - opt$recomb_rate*0.01/1e6)^(tmp_data$dist))/2
  # Bin data
  if(opt$fit_bin_size > 1) {
    breaks <- seq(0, max(tmp_data$dist)+opt$fit_bin_size, opt$fit_bin_size)
    tmp_data$dist <- cut(tmp_data$dist, breaks, head(breaks, -1))
    if(opt$bin_quant > 0) {
      tmp_data <- aggregate(. ~ dist, data=tmp_data, quantile, probs=opt$bin_quant/100, na.rm=TRUE)
    } else {
      tmp_data <- aggregate(. ~ dist, data=tmp_data, mean, na.rm=TRUE)
    }
  }
  tmp_data$File <- ld_file
  ld_data <- rbind(ld_data, melt(tmp_data, c("File","dist"), variable.name="LD", na.rm=TRUE))
}
# Clean-up
rm(tmp_data)
# Remove factor in dist
ld_data$dist <- as.numeric(levels(ld_data$dist))[ld_data$dist]
# Set maximum X-axis
if(is.null(opt$plot_x_lim)) {
  opt$plot_x_lim = max(ld_data$dist)
} else {
  opt$plot_x_lim = opt$plot_x_lim * 1000
}
# Set maximum Y-axis
if(!is.null(opt$plot_y_lim))
  opt$plot_y_lim <- c(0, opt$plot_y_lim)
# Add extra info
ld_data <- merge(ld_files, ld_data, by="File", sort=FALSE)
# 
n_groups <- length(unique(ld_data[,opt$plot_group]))
n_plots <- n_files * n_ld / n_groups
lgr$debug("Final data", data = toJSON(head(ld_data)))



### Fit decay
# Model function
ld_exp <- function(par, dist, ld_stat) {
  par <- as.numeric(par)
  if(ld_stat == "r2" || ld_stat == "r2_ExpG") {
    C = par[1] * dist
    r2h = par[2]
    r2l = par[3]
    if(opt$n_ind){
      # LD decay curve adjusted for finite samples sizes
      ((10+C) / ((2+C)*(11+C))) * (1+((3+C)*(12+12*C+C^2))/(opt$n_ind*(2+C)*(11+C)))
    }else{
      # Theoretical expectation under to drift
      #1 / (1 + C)
      # Theoretical expectation with r2_high and r2_low
      (r2h - r2l) / (1 + C) + r2l
    }
  } else if(ld_stat == "Dp") {
    D0 = 1
    t = par[1]
    Dh = par[2]
    Dl = par[3]
    Dl + (Dh-Dl) * D0 * (1 - dist * opt$recomb_rate/1e6)^t
  }
}
# Evaluation function
fit_eval <- function(par, obs_data) {
  if(length(unique(obs_data$LD)) != 1)
    stop(lgr$error("Invalid data.frame (several LD measures)"))
  model <- ld_exp(par, obs_data$dist, obs_data$LD[1])
  eval <- sum((model - obs_data$value)^2)
  lgr$trace("Fit evaluation", pars = par, eval = eval)
  eval
}
# Fitting function
fit_func <- function(x, fit_level) {
  ld_stat <- x$LD[1]
  # There is no fitting model for D
  if(ld_stat == "D") return(NULL)

  optim_tmp <- list()
  n_iter <- ifelse(fit_level>=10,fit_level,1)
  for(iter in 1:n_iter){
    # Fit LD model
    init_vals <- runif(3)
    if(ld_stat == 'Dp') {
      init_vals[1] = runif(1,10,20)
      par <- data.frame(init=init_vals, low_lim=c(0,0,0), up_lim=c(Inf,1,1))
    } else { # r2 and r2_ExpG
      init_vals[1] = runif(1,0,0.1)
      par <- data.frame(init=init_vals, low_lim=c(0,0,0), up_lim=c(1,1,1))
    }
    lgr$trace("Starting parameters", par = toJSON(par))
    
    optim_tmp <- append(optim_tmp, list("BFGS" = optim(par$init, fit_eval, obs_data=x, method="BFGS")) )
    if(fit_level > 1) optim_tmp <- append(optim_tmp, list("Nelder-Mead" = optim(par$init, fit_eval, obs_data=x, method="Nelder-Mead")) )
    if(fit_level > 2) optim_tmp <- append(optim_tmp, list("L-BFGS-B" = optim(par$init, fit_eval, obs_data=x, method="L-BFGS-B", lower=par$low_lim, upper=par$up_lim)) )
  }
  lgr$debug("Optimization", data = toJSON(optim_tmp))

  # If not using the theoretical r2 decay curve (with r2h and r2l)
  if(opt$n_ind > 0)
    if(ld_stat != 'Dp') optim_tmp <- lapply(optim_tmp, function(x){x$par[2]=x$par[3]=0;x})
  
  # Filter out runs that not-converged and/or with out-of-bound parameters
  # Columns stand for: par1 (rate), par2 (high), par3 (low), score, counts, convergence, message
  optim_tmp <- Filter(function(x){x$convergence == 0 &
                                  x$par[1] >= par$low_lim[1] & x$par[1] <= par$up_lim[1] & 
                                  x$par[2] >= par$low_lim[2] & x$par[2] <= par$up_lim[2] &
                                  x$par[3] >= par$low_lim[3] & x$par[3] <= par$up_lim[3] &
                                  x$par[2] >= x$par[3]}, optim_tmp)
  
  # Pick best run
  optim_tmp <- optim_tmp[order(sapply(optim_tmp,'[[',2))[1]]
  if(length(optim_tmp[[1]]) == 0) stop(lgr$error("convergence analyses failed. Please try increasing the fit level (--fit_level)"))
  lgr$debug("Best fit for %s (%s): %s", as.character(x$File[1]), as.vector(ld_stat), names(optim_tmp))
  optim_tmp[[1]]$par
}

# Fit LD decay distribution
if(opt$fit_level > 0) {
  # Define line "resolution"
  smooth <- seq(1, opt$plot_x_lim, length=opt$plot_line_smooth)
  # Full data
  optim_fit <- ddply(ld_data, .(File,LD), fit_func, fit_level=opt$fit_level)

  # Bootstrap
  if(opt$fit_boot > 0) {
    boot_rep_fit <- c()
    for (b in 1:opt$fit_boot) {
      ld_bootdata <- ddply(ld_data, .(File, LD), function(x) {x[sample(nrow(x), size=nrow(x), replace=TRUE),]} )
      boot_rep_fit <- rbind(boot_rep_fit, data.frame(ddply(ld_bootdata, .(File,LD), fit_func, fit_level=opt$fit_level), Rep=b))
    }
    boot_fit <- as.data.frame(as.matrix(aggregate(cbind(V1,V2,V3) ~ File+LD, boot_rep_fit, quantile, probs=c(0.025,0.975), names=FALSE)), stringsAsFactors=FALSE)
    optim_fit <- merge(optim_fit, boot_fit, sort=FALSE)
    optim_data <- ddply(optim_fit, .(File,LD), function(x) data.frame(LD=x[,"LD"], dist=smooth, value=ld_exp(x[,c("V1","V2","V3")], smooth, x[,"LD"]), ci_l=ld_exp(x[,c("V1.2","V2.1","V3.1")], smooth, x[,"LD"]), ci_u=ld_exp(x[,c("V1.1","V2.2","V3.2")], smooth, x[,"LD"])) )
  } else {
    optim_data <- ddply(optim_fit, .(File,LD), function(x) data.frame(LD=x[,"LD"], dist=smooth, value=ld_exp(x[,c("V1","V2","V3")], smooth, x[,"LD"])) )
  }
  # Print best FIT parameters
  optim_fit <- plyr::rename(optim_fit, c("V1"="DecayRate","V2"="LDmax","V3"="LDmin","V1.1"="DecayRate_CI.u","V1.2"="DecayRate_CI.l","V2.1"="LDmax_CI.l","V2.2"="LDmax_CI.u","V3.1"="LDmin_CI.l","V3.2"="LDmin_CI.u"), warn_missing=FALSE)
  lgr$info("Optimized fit", data = toJSON(optim_fit))

  # Merge data together with extra info from input
  fit_data <- merge(ld_files, optim_data, sort=FALSE)
  lgr$debug("Fit data", data = toJSON(head(fit_data, n=10)))
}



### Create base plot
lgr$info("Plotting data...")
plot <- ggplot() + 
  theme(panel.spacing=unit(1,"lines")) +
  coord_cartesian(xlim=c(0,opt$plot_x_lim), ylim=opt$plot_y_lim) +
  scale_colour_hue() +
  ylab("Linkage Disequilibrium") +
  xlab("Distance")

if(!is.null(opt$plot_wrap_formula)) {
  if(opt$plot_wrap) {
    plot <- plot + facet_wrap(opt$plot_wrap_formula, ncol=opt$plot_wrap, scales=opt$plot_axis_scales)
  } else {
    plot <- plot + facet_grid(opt$plot_wrap_formula, scales=opt$plot_axis_scales)
  }
}

# Add LD decay fit CI
if(opt$fit_boot > 0) {
  grp <- NULL
  if(n_files == n_groups)
    grp <- opt$plot_group
  plot <- plot + geom_ribbon(data=fit_data, aes(x=dist,ymin=ci_l,ymax=ci_u,group=.data[[opt$plot_group]],fill=.data[[grp]]), alpha=0.2)
}

# Add data points
if(opt$plot_data){
  # Check format
  if(ncol(ld_data) < 4)
    stop(lgr$error("Invalid `ld_data` format."))
  # Bin data
  if(opt$plot_bin_size > 1) {
    breaks <- seq(0, max(ld_data$dist)+opt$plot_bin_size, opt$plot_bin_size)
    ld_data$dist <- cut(ld_data$dist, breaks, head(breaks, -1))
    ld_data$dist <- as.numeric(levels(ld_data$dist))[ld_data$dist]
    if(opt$bin_quant > 0) {
      ld_data <- aggregate(value ~ ., data=ld_data, quantile, probs=opt$bin_quant/100)
    } else {
      ld_data <- aggregate(value ~ ., data=ld_data, mean)
    }
  }
  lgr$debug("Plot data", data = toJSON(head(ld_data, n=10)))
  # Add points
  plot <- plot + geom_point(data=ld_data, aes(x=dist,y=value,colour=.data[[opt$plot_group]]), size=0.05, alpha=0.2)
}

# Add LD decay best fit
if(length(opt$ld) > 0) {
  # Select fields that are variable
  header <- names(which(lapply(lapply(fit_data, unique), length) > 1))
  # Exclude non-relevant fields
  grp <- header[!header %in% unique(c(as.character(opt$plot_wrap_formula), opt$plot_group, "dist", "value", "File", "ci_l", "ci_u"))]
  lgr$debug("LD decay best fit", groups = grp, plot_groups = opt$plot_group)
  # Define line type
  if(length(grp) == 0) grp <- opt$plot_group
  if(length(grp) > 1) stop(lgr$error("invalid number of linetype groups!"))
  plot <- plot + geom_line(data=fit_data, aes(x=dist,y=value,colour=.data[[opt$plot_group]],linetype=.data[[grp]]))
  # If ploting data, add a thin black line to help see the line
  if(opt$plot_data && FALSE) # Disabled since, when plottig more than 1 line, it plots a grey area between them
    plot <- plot + geom_line(data=fit_data, aes(x=dist,y=value,colour=.data[[opt$plot_group]]), size=0.1, alpha=0.2, colour="black")
  # If plotting, apart from linetypes, also shapes (for B/W or color-blind printing)
  if(opt$plot_shape) {
    smooth <- seq(1, opt$plot_x_lim, length=opt$plot_line_smooth)[seq(2,opt$plot_line_smooth,length=5)]
    sample_fit_data <- subset(fit_data, dist %in% smooth)
    plot <- plot + geom_point(data=sample_fit_data, aes(x=dist,y=value,colour=.data[[opt$plot_group]],shape=.data[[opt$plot_group]]))
  }
}



### Set plot size
n_plots <- length(unique(ggplot_build(plot)$data[[1]]$PANEL))
if(!is.null(opt$plot_wrap_formula)) {
  par <- dcast(fit_data, opt$plot_wrap_formula, length, fill=0)
  rownames(par) <- par[,1]
  par <- par[,-1, drop=FALSE]
} else {
  par <- matrix(ncol=1)
}
lgr$debug("Summary", n_files = n_files, n_ld = n_ld, n_groups = n_groups, n_plots = n_plots)
lgr$debug("Parameters layout", rows_par = nrow(par), cols_par = ncol(par))
plot_height <- opt$plot_size[1] * nrow(par)
plot_width <- opt$plot_size[2] * ncol(par)



### Remove legend if plotting just a single variable
if(n_groups < 2 || opt$plot_no_legend) {
  plot <- plot + theme(legend.position="none")
} else {
  plot_width = plot_width + 1
}



### Save plot
ggsave(opt$out, plot=plot, height=plot_height, width=plot_width, scale=opt$plot_scale, limitsize=FALSE)
