#!/bin/env Rscript

#FileName: fit_LDdecay.R v1.0
#Author: "Filipe G. Vieira (fgarrettvieira _at_ gmail [dot] com)"
#Author: "Emma Fox (e.fox16 _at_ imperial [dot] ac [dot] uk)"

##IMPORTS
library(optparse)
library(tools)
library(ggplot2)
library(reshape2)
library(plyr)
options(width = 160) 



##SPECIFY ARGUMENTS
option_list <- list(
  make_option(c('--ld_files'), action='store', type='character', default=NULL, help = 'File with list of LD files to fit and plot (if ommited, can be read from STDIN)'),
  make_option(c('--header'), action='store_true', type='logical', default=FALSE, help='Input file has header'),
  make_option(c('--col'), action='store', type='numeric', default=3, help='Which column is distance between sites? [%default]'),
  make_option(c('--ld'), action='store', type='character', default="r2", help='Which LD stats to plot (r2pear, D, Dp, r2) [%default]'),
  make_option(c('--n_ind'), action='store', type='numeric', default=0, help='How many individuals in the sample (for r^2 fitting correction)?'),
  make_option(c('--max_kb_dist'), action='store', type='numeric', default=Inf, help='Maximum distance between SNPs (in kb) to include in the fitting analysis. [%default]'),
  make_option(c('--fit_boot'), action='store', type='numeric', default=0, help='Number of bootstrap replicates for fitting CI. [%default]'),
  make_option(c('--fit_bin_size'), action='store', type='numeric', default=250, help='Bin data into fixed-sized windows and use the average for fitting. [default %default bps]'),
  make_option(c('--fit_level'), action='store', type='numeric', default=1, help='Fitting level 0) no fitting, best of 1) Nelder-Mead, 2) and BFGS, 3) and L-BFGS-B). [%default]'),
  make_option(c('--plot_group'), action='store', type='character', default='File', help='Group variable'),
  make_option(c('--plot_data'), action='store_true', type='logical', default=FALSE, help='Also plot data points?'),
  make_option(c('--plot_bin_size'), action='store', type='numeric', default=1000, help='Bin data into fixed-sized windows and use the average for plotting. [default %default bps]'),
  make_option(c('--plot_x_lim'), action='store', type='numeric', default=0, help='X-axis plot limit (in kb). [%default]'),
  make_option(c('--plot_axis_scales'), action='store', type='character', default='fixed', help='Plot axis scales: fixed (default), free, free_x or free_y'),
  make_option(c('--plot_size'), action='store', type='character', default='1,1', help='Plot size (height,width). [%default]'),
  make_option(c('--plot_scale'), action='store', type='numeric', default=2.5, help='Plot scale. [%default]'),
  make_option(c('--plot_wrap'), action='store', type='numeric', default=0, help='Plot in WRAP with X columns (default in GRID)'),
  make_option(c('--plot_no_legend'), action='store_true', type='logical', default=FALSE, help='REmove legend from plot'),
  make_option(c('-f','--plot_wrap_formula'), action='store', type='character', default=NULL, help='Plot formula for WRAP. [%default]'),
  make_option(c('-o','--out'), action='store', type='character', default=NULL, help='Output file'),
  make_option(c('--debug'), action='store_true', type='logical', default=FALSE, help='Debug mode. Extra output...')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Parse input LD files
if(is.null(opt$ld_files)) 
  opt$ld_files <- file('stdin')
ld_files <- read.table(opt$ld_files, header=opt$header, stringsAsFactors=FALSE)
colnames(ld_files)[1] <- "File"
# Keep file list ordered
for (id in colnames(ld_files))
    ld_files[,id] <- factor(ld_files[,id], unique(ld_files[,id]), ordered=TRUE)
n_files <- nrow(ld_files)
if(opt$debug)
  print(ld_files)

# Parse LD stats to plot
if(!is.null(opt$ld))
  opt$ld <- unlist(strsplit(opt$ld, ","))
if(!all(opt$ld %in% c("r2pear", "D", "Dp", "r2")))
  stop("Invalid LD measure to plot")
n_ld = length(opt$ld)

# Check if number of individuals was specified
if(any(opt$ld %in% c("r2pear", "r2")) && !opt$n_ind)
  stop("Fitting of R^2 requires number of individuals")

# Check max_kb_dist parameter 
if(opt$max_kb_dist < 50)
  warning("Fitting of LD decay is highly unreliable at short distances (<50kb).")

# Check binning sizes
if(opt$fit_bin_size > opt$plot_bin_size)
   stop("Fitting bin size cannot be greater than plotting bin size!")

# Parse plot_size
opt$plot_size <- as.numeric(unlist(strsplit(opt$plot_size, ",")))
if(length(opt$plot_size) < 2)
  opt$plot_size <- c(opt$plot_size, opt$plot_size)

# Plot formula
if(!is.null(opt$plot_wrap_formula))
  opt$plot_wrap_formula <- as.formula(opt$plot_wrap_formula)

# Set output file name (if not defined)
if(is.null(opt$out))
  opt$out <- paste(basename(file_path_sans_ext(opt$ld_files)),".pdf", sep = "")



### Load LD data
header <- c("Dist","r2pear","D","Dp","r2")
ld_data <- data.frame()
for (i in 1:n_files) {
  ld_file <- as.character(ld_files$File[i])
  # Read point data
  if(opt$debug)
    cat("Reading file:", ld_file, fill=TRUE)
  tmp_data <- read.table(gzfile(ld_file), sep="\t", quote="\"", dec=".")[-(1:(opt$col-1))]
  if(ncol(tmp_data) < 5)
    stop('Invalid LD file format.\n')
  colnames(tmp_data) <- header
  tmp_data <- tmp_data[, which(names(tmp_data) %in% c("Dist",opt$ld))]
  # Filter by minimum distance
  tmp_data <- tmp_data[which(tmp_data$Dist < opt$max_kb_dist*1000),]
  # Bin data
  if(opt$fit_bin_size != 0) {
    tmp_data$Dist <- as.integer(tmp_data$Dist / opt$fit_bin_size) * opt$fit_bin_size
    tmp_data <- aggregate(. ~ Dist, data=tmp_data, median)
  }
  tmp_data$File <- ld_file
  ld_data <- rbind(ld_data, melt(tmp_data, c("File","Dist"), variable.name="LD", na.rm=TRUE))
}
# Clean-up
rm(tmp_data)
# Set maximum X-axis
opt$plot_x_lim = opt$plot_x_lim * 1000
if(opt$plot_x_lim == 0)
  opt$plot_x_lim = max(ld_data$Dist)
# Add extra info
ld_data <- merge(ld_files, ld_data, by="File", sort=FALSE)
# 
n_groups <- length(unique(ld_data[,opt$plot_group]))
n_plots <- n_files * n_ld / n_groups
if(opt$debug)
  print(head(ld_data))



### Fit decay
# Model function
ld_exp <- function(par, dist, ld_stat) {
  par <- as.numeric(par)
  if(ld_stat == "r2" || ld_stat == "r2pear") {
    C = par[1] * dist
    r2h = par[2]
    r2l = par[3]
    if(0){
      # Theoretical expectation
      #1 / (1 + C)
      # Theoretical expectation with Dh and Dl (not recomended since it assumes an infinite sample size)
      (r2h - r2l) / (1 + C) + r2l
    }else{
      # Adjusted for finite samples sizes
      ((10+C) / ((2+C)*(11+C))) * (1+((3+C)*(12+12*C+C^2))/(opt$n_ind*(2+C)*(11+C)))
    }
  } else if(ld_stat == "Dp") {
    D0 = 1
    t = par[1]
    Dh = par[2]
    Dl = par[3]
    Dl + (Dh-Dl) * D0 * (1 - dist/1000000)^t # Assuming 1cm = 1Mb
  }
}
# Evaluation function
fit_eval <- function(par, obs_data) {
  if(length(unique(obs_data$LD)) != 1)
    stop("Invalid data.frame (several LD measures)")
  model <- ld_exp(par, obs_data$Dist, obs_data$LD[1])
  sum((model - obs_data$value)^2)
}
# Fitting function
fit_func <- function(x, fit_level) {
  ld_stat <- x$LD[1]
  # There is no fitting model for D
  if(ld_stat == "D") return(NULL)
  
  # Fit LD model
  if(ld_stat == 'Dp') {
    par <- data.frame(init=c(10,0.9,0.1), low_lim=c(0,0,0), up_lim=c(Inf,1,1))
  } else { # r2 and r2pear
    par <- data.frame(init=c(0.01,0.9,0.1), low_lim=c(0,0,0), up_lim=c(1,1,1))
  }
  optim_tmp <- list("BFGS" = optim(par$init, fit_eval, obs_data=x, method="BFGS"))
  if(fit_level > 1) optim_tmp <- append(optim_tmp, list("Nelder-Mead" = optim(par$init, fit_eval, obs_data=x, method="Nelder-Mead")) )
  if(fit_level > 2) optim_tmp <- append(optim_tmp, list("L-BFGS-B" = optim(par$init, fit_eval, obs_data=x, method="L-BFGS-B", lower=par$low_lim, upper=par$up_lim)) )
  if(opt$debug) lapply(optim_tmp, function(x){ cat(format(as.numeric(x$par),digits=4,nsmall=6,scientific=FALSE), format(as.numeric(x$value),digits=7,nsmall=4), x$counts, x$convergence, x$message, sep="\t", fill=TRUE) })
  # If not using the theoretical D' decay curve (with Dh and Dl)
  if(!0)
    if(ld_stat != 'Dp') optim_tmp <- lapply(optim_tmp, function(x){x$par[2]=x$par[3]=0;x})
  
  # Filter out runs that not-converged and/or with out-of-bound parameters
  # Columns stand for: par1, par2, par3, score, counts, convergence, message
  optim_tmp <- Filter(function(x){x$convergence == 0 &
                                  x$par[1] >= par$low_lim[1] & x$par[1] <= par$up_lim[1] & 
                                  x$par[2] >= par$low_lim[2] & x$par[2] <= par$up_lim[2] &
                                  x$par[3] >= par$low_lim[3] & x$par[3] <= par$up_lim[3]}, optim_tmp)
  
  # Pick best run
  optim_tmp <- optim_tmp[order(sapply(optim_tmp,'[[',2))[1]]
  if(length(optim_tmp) != 1) stop("convergence analyses failed.")
  if(opt$debug) cat("Best fit for ", as.character(x$File[1]), " (",as.vector(ld_stat),"): ", names(optim_tmp), sep="", fill=TRUE)
  optim_tmp[[1]]$par
}

# Fit LD decay distribution
if(opt$fit_level > 0) {
  grid <- seq(1, opt$plot_x_lim, length=1000)
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
    all_fit <- merge(optim_fit, boot_fit, sort=FALSE)
    optim_data <- ddply(all_fit, .(File,LD), function(x) data.frame(LD=x[,"LD"], Dist=grid, value=ld_exp(x[,c("V1","V2","V3")], grid, x[,"LD"]), ci_l=ld_exp(x[,c("V1.2","V2.1","V3.1")], grid, x[,"LD"]), ci_u=ld_exp(x[,c("V1.1","V2.2","V3.2")], grid, x[,"LD"])) )
  } else {
    optim_data <- ddply(optim_fit, .(File,LD), function(x) data.frame(LD=x[,"LD"], Dist=grid, value=ld_exp(x[,c("V1","V2","V3")], grid, x[,"LD"])) )
  }
  # Merge data together with extra info from input
  fit_data <- merge(ld_files, optim_data, sort=FALSE)
}
# Print 
print(optim_fit)
if(opt$debug)
  print(head(fit_data))



### Create base plot
plot <- ggplot() + xlim(0, opt$plot_x_lim) + theme(panel.spacing=unit(1,"lines"))
if(!is.null(opt$plot_wrap_formula)) {
  if(opt$plot_wrap) {
    #cat('# WRAP mode', fill=TRUE)
    plot <- plot + facet_wrap(opt$plot_wrap_formula, ncol=opt$plot_wrap, scales=opt$plot_axis_scales)
  } else {
    #cat('# GRID mode', fill=TRUE)
    plot <- plot + facet_grid(opt$plot_wrap_formula, scales=opt$plot_axis_scales)
  }
}

# Add LD decay fit CI
if(opt$fit_boot > 0)
  plot <- plot + geom_ribbon(data=fit_data, aes(x=Dist,ymin=ci_l,ymax=ci_u), alpha=0.2)

# Add data points
if(opt$plot_data){
  # Check format
  if(ncol(ld_data) < 4)
    stop("Invalid `ld_data` format.\n")
  # Bin data
  if(opt$plot_bin_size > 1) {
    ld_data$Dist <- as.integer(ld_data$Dist / opt$plot_bin_size) * opt$plot_bin_size
    ld_data <- aggregate(value ~ ., data=ld_data, median)
  }
  # Add points
  plot <- plot + geom_point(data=ld_data, aes_string(x="Dist",y="value",colour=opt$plot_group), size=0.05, alpha=0.2)
}

# Add LD decay best fit
if(length(opt$ld) > 0) {
  header <- colnames(fit_data)
  grp <- header[!header %in% unique(c(as.character(opt$plot_wrap_formula), opt$plot_group, "Dist", "value", "File", "ci_l", "ci_u"))]
  if(length(grp) == 0) grp <- NULL
  plot <- plot + geom_line(data=fit_data, aes_string(x="Dist",y="value",group=opt$plot_group,colour=opt$plot_group,linetype=grp))
  if(opt$plot_data)
    plot <- plot + geom_line(data=fit_data, aes_string(x="Dist",y="value",group=opt$plot_group), size=0.1, alpha=0.2, colour="black")
}



### Set plot size
n_plots <- length(unique(ggplot_build(plot)$data[[1]]$PANEL))
if(!is.null(opt$plot_wrap_formula)) {
  par <- dcast(fit_data, opt$plot_wrap_formula, length, fill=0)
  rownames(par) <- par[,1]; par <- par[,-1]
} else {
  par <- matrix()
}
if(opt$debug){
  cat(n_files, n_ld, n_groups, n_plots, fill=TRUE)
  cat(nrow(par), ncol(par), fill=TRUE)
}
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
x <- warnings()
