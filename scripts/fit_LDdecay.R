#!/usr/bin/Rscript

#FileName: fit_LDdecay.R
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
  make_option(c('--col'), action='store', type='numeric', default=3, help='Which column is distance between sites'),
  make_option(c('--ld'), action='store', type='character', default="r2", help='Which LD stats to plot (r2pear, D, Dp, r2)'),
  make_option(c('--max_dist'), action='store', type='numeric', default=Inf, help='Maximum distance between SNPs (in kb) to include in the fitting analysis'),
  make_option(c('--fit_boot'), action='store', type='numeric', default=0, help='Number of bootstrap replicates for fitting CI'),
  make_option(c('--fit_bin_size'), action='store', type='numeric', default=0, help='Size of bin to fit'),
  make_option(c('--fit_level'), action='store', type='numeric', default=1, help='Fitting level 0) no fitting, best of 1) Nelder-Mead, 2) and BFGS, 3) and L-BFGS-B)'),
  make_option(c('--plot_group'), action='store', type='character', default='File', help='Group variable'),
  make_option(c('--plot_data'), action='store_true', type='logical', default=FALSE, help='Also plot data points?'),
  make_option(c('--plot_bin_size'), action='store', type='numeric', default=100, help='Size of bin to plot'),
  make_option(c('--plot_x_lim'), action='store', type='numeric', default=1e6, help='X-axis plot limit'),
  make_option(c('--plot_axis_scales'), action='store', type='character', default='fixed', help='Plot axis scales: fixed (default), free, free_x or free_y'),
  make_option(c('--plot_size'), action='store', type='character', default='1,1', help='Plot size (height,width)'),
  make_option(c('--plot_scale'), action='store', type='numeric', default=2, help='Plot scale'),
  make_option(c('--plot_wrap'), action='store', type='numeric', default=0, help='Plot in WRAP with X columns (default in GRID)'),
  make_option(c('-f','--plot_wrap_formula'), action='store', type='character', default='File ~ LD', help='Plot formula for WRAP'),
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
n_files <- nrow(ld_files)
if(opt$debug)
  print(ld_files)

# Parse LD stats to plot
if(!is.null(opt$ld))
  opt$ld <- unlist(strsplit(opt$ld, ","))
if(!all(opt$ld %in% c("r2pear", "D", "Dp", "r2")))
  error("Invalid LD measure to plot")
n_ld = length(opt$ld)

# Check max_dist parameter 
if(opt$max_dist < 50)
  warning("Fitting of LD decay is highly unreliable at short distances (<50kb).")

# If there is binning for the fit, plot that data (do not bin again)
if(opt$fit_bin_size)
  opt$plot_bin_size = 0

# Parse plot_size
opt$plot_size <- as.numeric(unlist(strsplit(opt$plot_size, ",")))
if(length(opt$plot_size) < 2)
  opt$plot_size <- c(opt$plot_size, opt$plot_size)

# Plot formula
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
  tmp_data <- tmp_data[which(tmp_data$Dist < opt$max_dist*1000),]
  # Bin data
  if(opt$fit_bin_size != 0) {
    tmp_data$Dist <- as.integer(tmp_data$Dist / opt$fit_bin_size) * opt$fit_bin_size + 1
    tmp_data <- aggregate(. ~ Dist, data=tmp_data, mean)
  }
  tmp_data$File <- ld_file
  ld_data <- rbind(ld_data, melt(tmp_data, c("File","Dist"), variable.name="LD", na.rm=TRUE))
}
# Clean-up
rm(tmp_data)
# Add extra info
ld_data <- merge(ld_data, ld_files, by="File")
# 
n_groups <- length(unique(ld_data[,opt$plot_group]))
n_plots <- n_files / n_groups
if(opt$debug) {
  cat(n_files, n_groups, n_ld, n_plots, fill=TRUE)
  print(head(ld_data))
}


### Fit decay
# Model function
ld_exp <- function(par, dist, ld_stat) {
  par <- as.numeric(par)
  if(ld_stat == "r2" || ld_stat == "r2pear") {
    #1 / (1 + par[1] * dist) # Theoretical expectation
    C = par[1] * dist
    n = 50
    ((10+C) / ((2+C)*(11+C))) * (1+((3+C)*(12+12*C+C^2))/(n*(2+C)*(11+C)))
  } else if(ld_stat == "Dp") {
    D0 = 1
    Dl = par[1]
    Dh = par[2]
    t = par[3]
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
    par <- data.frame(init=c(0.01,0.9,10), low_lim=c(0,0,0), up_lim=c(1,1,Inf))
    optim_tmp <- list("BFGS" = optim(par$init, fit_eval, obs_data=x, method=c("BFGS")))
    if(fit_level > 1) optim_tmp <- append(optim_tmp, list("Nelder-Mead" = optim(par$init, fit_eval, obs_data=x, method=c("Nelder-Mead"))) )
    if(fit_level > 2) optim_tmp <- append(optim_tmp, list("L-BFGS-B" = optim(par$init, fit_eval, obs_data=x, method=c("L-BFGS-B"), lower=par$low_lim, upper=par$up_lim)) )
  } else { # r2 and r2pear
    par <- data.frame(init=c(0.01), low_lim=c(0), up_lim=c(1))
    optim_tmp <- list("Brent" = optim(par$init, fit_eval, obs_data=x, method=c("Brent"), lower=par$low_lim, upper=par$up_lim))
    optim_tmp[[1]]$par <- c(optim_tmp[[1]]$par, 0, 0)
  }
  if(opt$debug) lapply(optim_tmp, function(x){ cat(format(as.numeric(x$par),digits=4,nsmall=6,scientific=FALSE), format(as.numeric(x$value),digits=7,nsmall=4), x$counts, x$convergence, x$message, sep="\t", fill=TRUE) } )
  # Filter out runs that not-converged and/or with out-of-bound parameters
  optim_tmp <- Filter(function(x){x$convergence == 0 &
                                  x$par[1] >= par$low_lim[1] & x$par[1] <= par$up_lim[1] & 
                                  (ld_stat != 'Dp' | x$par[2] >= par$low_lim[2] & x$par[2] <= par$up_lim[2] &
                                                     x$par[3] >= par$low_lim[3] & x$par[3] <= par$up_lim[3])}, optim_tmp)
  # Pick best run
  optim_tmp <- optim_tmp[order(sapply(optim_tmp,'[[',2))[1]]
  if(length(optim_tmp) != 1) stop("convergence analyses failed.")
  if(opt$debug) cat("Best fit for ", x$File[1], " (",as.vector(ld_stat),"): ", names(optim_tmp), sep="", fill=TRUE)
  optim_tmp[[1]]$par
}

# Fit LD decay distribution
if(opt$fit_level > 0) {
  grid <- seq(1, opt$plot_x_lim, length=1000)
  # Full data
  optim_fit <- ddply(ld_data, .(File,LD), fit_func, fit_level=opt$fit_level)
  if(opt$debug)
    print(optim_fit)
  # Bootstrap
  if(opt$fit_boot > 0) {
    boot_rep_fit <- c()
    for (b in 1:opt$fit_boot) {
      ld_bootdata <- ddply(ld_data, .(File, LD), function(x) {x[sample(nrow(x), size=nrow(x), replace=TRUE),]} )
      boot_rep_fit <- rbind(boot_rep_fit, data.frame(ddply(ld_bootdata, .(File,LD), fit_func, fit_level=opt$fit_level), Rep=b))
    }
    boot_fit <- as.data.frame(as.matrix(aggregate(cbind(V1,V2,V3) ~ File+LD, boot_rep_fit, quantile, probs=c(0.025,0.975), names=FALSE)), stringsAsFactors=FALSE)
    all_fit <- merge(optim_fit, boot_fit)
    optim_data <- ddply(all_fit, .(File,LD), function(x) data.frame(LD=x[,"LD"], Dist=grid, value=ld_exp(x[,c("V1","V2","V3")], grid, x[,"LD"]), ci1=ld_exp(x[,c("V1.1","V2.1","V3.1")], grid, x[,"LD"]), ci2=ld_exp(x[,c("V1.2","V2.2","V3.2")], grid, x[,"LD"])) )
  } else {
    optim_data <- ddply(optim_fit, .(File,LD), function(x) data.frame(LD=x[,"LD"], Dist=grid, value=ld_exp(x[,c("V1","V2","V3")], grid, x[,"LD"])) )
  }
  # Merge data together with exrta info from input
  fit_data <- merge(optim_data, ld_files)
}
# Print 
print(optim_fit)
if(opt$debug)
  print(head(fit_data))



### Create base plot
plot <- ggplot() + xlim(0, opt$plot_x_lim)
if(opt$plot_wrap) {
  #cat('# WRAP mode', fill=TRUE)
  plot <- plot + facet_wrap(opt$plot_wrap_formula, ncol=opt$plot_wrap, scales=opt$plot_axis_scales)
} else {
  #cat('# GRID mode', fill=TRUE)
  plot <- plot + facet_grid(opt$plot_wrap_formula, scales=opt$plot_axis_scales)
}
# Add LD decay fit CI
if(opt$fit_boot > 0)
  plot <- plot + geom_ribbon(data=fit_data, aes(x=Dist,ymin=ci1,ymax=ci2), alpha=0.2)
# Add data points
if(opt$plot_data){
  # Check format
  if(ncol(ld_data) < 4)
    stop("Invalid `ld_data` format.\n")
  # Bin data
  if(opt$plot_bin_size != 0) {
    ld_data$Dist <- as.integer(ld_data$Dist / opt$plot_bin_size) * opt$plot_bin_size
    ld_data <- aggregate(value ~ ., data=ld_data, mean)
  }
  # Add points
  plot <- plot + geom_point(data=ld_data, aes_string(x="Dist",y="value",colour=opt$plot_group), size=0.1, alpha=0.2)
}
# Add LD decay best fit
if(length(opt$ld) > 0) {
  header <- colnames(fit_data)
  grp <- header[!header %in% unique(c(as.character(opt$plot_wrap_formula), opt$plot_group, "Dist", "value", "File", "ci1", "ci2"))]
  if(length(grp) == 0) grp <- NULL
  plot <- plot + geom_line(data=fit_data, aes_string(x="Dist",y="value",group=opt$plot_group,colour=opt$plot_group,linetype=grp))
}



### Save plot
plot_height <- opt$plot_size[1] * n_plots
plot_width <- opt$plot_size[2] * n_groups
ggsave(opt$out, plot=plot, height=plot_height, width=plot_width, scale=opt$plot_scale, limitsize=FALSE)
x <- warnings()
