#!/usr/bin/Rscript

#FileName: fit_LDdecay.R
#Author: "Filipe G. Vieira (fgarrettvieira _at_ gmail [dot] com)"
#Author: "Emma Fox (e.fox16 _at_ imperial [dot] ac [dot] uk)"

##IMPORTS
library(optparse)
library(tools)
library(ggplot2)
library(reshape2)
library(dplyr)
#library(optimx)



##SPECIFY ARGUMENTS
option_list <- list(
  make_option(c('--ld_files'), action='store', type='character', default=NULL, help = 'File with list of LD files to fit and plot (if ommited, can be read from STDIN)'),
  make_option(c('--header'), action='store_true', type='logical', default=FALSE, help='Input file has header'),
  make_option(c('--col'), action='store', type='numeric', default=3, help='Which column is distance between sites'),
  make_option(c('--ld'), action='store', type='character', default="r2", help='Which LD stats to plot (r2pear, D, Dp, r2)'),
  make_option(c('--max_dist'), action='store', type='numeric', default=Inf, help="Maximum distance between SNPs to include in the fitting analysis"),
  make_option(c('--fit_bin_size'), action='store', type='numeric', default=0, help='Size of bin to fit'),
  make_option(c('--plot_group'), action='store', type='character', default='File', help='Group variable'),
  make_option(c('--plot_data'), action='store_true', type='logical', default=FALSE, help='Also plot data points?'),
  make_option(c('--plot_bin_size'), action='store', type='numeric', default=100, help='Size of bin to plot'),
  make_option(c('--plot_x_lim'), action='store', type='numeric', default=1e6, help='X-axis plot limit'),
  make_option(c('--plot_axis_scales'), action='store', type='character', default='fixed', help='Plot axis scales: fixed (default), free, free_x or free_y'),
  make_option(c('--plot_size'), action='store', type='character', default='1,2', help='Plot size (height,width)'),
  make_option(c('--plot_scale'), action='store', type='numeric', default=1, help='Plot scale'),
  make_option(c('--plot_wrap'), action='store', type='numeric', default=0, help='Plot in WRAP with X columns (default in GRID)'),
  make_option(c('-f','--plot_wrap_formula'), action='store', type='character', default='File ~ LD', help='Plot formula for WRAP'),
  make_option(c('-o','--out'), action='store', type='character', default=NULL, help='Output file'),
  make_option(c('--debug'), action='store_true', type='logical', default=FALSE, help='Debug mode. Extra output...')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#opt$ld_files = "Desktop/TEST/in.tsv"
#opt$ld = c("r2","Dp","r2pear","D")
#opt$header=TRUE
#opt$col = 5
#opt$plot_data = TRUE
#opt$debug = TRUE

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
n_ld = length(opt$ld)

# Check bin size options
if(opt$plot_bin_size && opt$fit_bin_size)
  stop('Cannot bin data twice.')
  
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
  tmp_data <- tmp_data[which(tmp_data$Dist < opt$max_dist),]
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
fit_data <- data.frame()
if(length(opt$ld) > 0) {
  # Model functions
  ld_exp <- function(par, dist, ld_stat) {
    if(ld_stat == "r2" || ld_stat == "r2pear") {
      par[1] / (1 + par[2] * dist) 
    } else if(ld_stat == "Dp") {
      stop("Dp fitting function not implemented yet.")
    } else {
      stop("Invalid LD stat specified")
    }
  }
  # Evaluation function
  fit_eval <- function(par, data, ld_stat) {
    data <- data[data$LD == ld_stat,]
    model <- ld_exp(par, data$Dist, ld_stat)
    sum((model - data$value)^2)
  }
  fit_func <- function(x, ld_stat) optim(c(0.1,0.0000001), fit_eval, data=ld_data[ld_data$File==x,], ld_stat=ld_stat, method="L-BFGS-B", lower=c(0,0), upper=c(1,1))
  #fit_func <- function(x, ld_stat) optimx(c(0.001,0.001), fit_eval, data=ld_data[ld_data$File==x,], ld_stat=ld_stat, method="BFGS")#, lower=c(0,0), upper=c(1,1))
  # Get fitted values
  grid <- seq(1, opt$plot_x_lim, length = 1000)
  for (ld_stat in opt$ld) {
    optim_fit <- as.data.frame(sapply(ld_files$File, fit_func, ld=ld_stat))
    optim_data <- lapply(optim_fit, function(x) data.frame(LD=ld_stat, Dist=grid, value=ld_exp(x$par, grid, ld_stat)) )
    fit_data <- rbind(fit_data, merge(bind_rows(optim_data, .id="File"), ld_files, by="File"))
    if(opt$debug)
      print(optim_fit)
  }
  if(opt$debug)
    print(head(fit_data))
}



### Create base plot
plot <- ggplot() + xlim(0, opt$plot_x_lim)
if(opt$plot_wrap) {
  #cat('# WRAP mode', fill=TRUE)
  plot <- plot + facet_wrap(opt$plot_wrap_formula, ncol=opt$plot_wrap, scales=opt$plot_axis_scales)
} else {
  #cat('# GRID mode', fill=TRUE)
  plot <- plot + facet_grid(opt$plot_wrap_formula, scales=opt$plot_axis_scales, space=opt$plot_axis_scales)
}
# Add data
if(opt$plot_data){
  # Check format
  if(ncol(ld_data) < 4)
    stop("Invalid `ld_data` format.\n")
  # Bin data
  if(opt$plot_bin_size != 0) {
    ld_data$Dist <- as.integer(ld_data$Dist / opt$plot_bin_size) * opt$plot_bin_size + 1
    ld_data <- aggregate(value ~ ., data=ld_data, mean)
  }
  # Add points
  plot <- plot + geom_point(data=ld_data, aes_string(x="Dist",y="value",colour=opt$plot_group), size=0.1, alpha=0.2)
}
# Add LD decay fit
if(length(opt$ld) > 0) {
  header <- colnames(fit_data)
  grp <-header[!header %in% unique(c(as.character(opt$plot_wrap_formula), opt$plot_group, "Dist", "value", "File"))]
  if(length(grp) == 0) grp <- NULL
  plot <- plot + geom_line(data=fit_data, aes_string(x="Dist",y="value",group=opt$plot_group,colour=opt$plot_group,linetype=grp))
}



### Save plot
plot_height <- 2 * opt$plot_size[1] * n_plots
plot_width <- 5 * opt$plot_size[2]
ggsave(opt$out, plot=plot, height=plot_height, width=plot_width, scale=opt$plot_scale, limitsize=FALSE)
