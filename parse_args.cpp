
#include <getopt.h>
#include "ngsLD.hpp"


void init_pars(params *pars) {
  pars->in_geno = NULL;
  pars->in_bin = false;
  pars->in_probs = false;
  pars->in_logscale = false;
  pars->n_ind = 0;
  pars->n_sites = 0;
  pars->in_pos = NULL;
  pars->in_pos_header = false;
  pars->max_kb_dist = 100;
  pars->max_snp_dist = 0;
  pars->min_maf = 0;
  pars->ignore_miss_data = false;
  pars->call_geno = false;
  pars->N_thresh = 0;
  pars->call_thresh = 0;
  pars->rnd_sample = 1;
  pars->seed = time(NULL) + rand() % 1000;
  pars->extend_out = false;
  pars->out = NULL;
  pars->out_fh = stdout;
  pars->n_threads = 1;
  pars->verbose = 1;
}


// Parses command line args and stores them into struct params
void parse_cmd_args(params* pars, int argc, char** argv) {
  
  static struct option long_options[] =
    {
      {"geno", required_argument, NULL, 'g'},
      {"probs", no_argument, NULL, 'p'},
      {"log_scale", no_argument, NULL, 'l'},
      {"n_ind", required_argument, NULL, 'n'},
      {"n_sites", required_argument, NULL, 's'},
      {"pos", required_argument, NULL, 'a'},
      {"posH", required_argument, NULL, 'A'},
      {"max_kb_dist", required_argument, NULL, 'd'},
      {"max_snp_dist", required_argument, NULL, 'D'},
      {"min_maf", required_argument, NULL, 'f'},
      {"ignore_miss_data", no_argument, NULL, 'm'},
      {"call_geno", no_argument, NULL, 'c'},
      {"N_thresh", required_argument, NULL, 'N'},
      {"call_thresh", required_argument, NULL, 'C'},
      {"rnd_sample", required_argument, NULL, 'r'},
      {"seed", required_argument, NULL, 'S'},
      {"extend_out", no_argument, NULL, 'x'},
      {"out", required_argument, NULL, 'o'},
      {"outH", required_argument, NULL, 'O'},
      {"n_threads", required_argument, NULL, 't'},
      {"verbose", required_argument, NULL, 'V'},
      {0, 0, 0, 0}
    };
  
  int c = 0;
  while ( (c = getopt_long_only(argc, argv, "g:pln:s:Z:d:D:f:mcN:C:r:S:xo:t:V:", long_options, NULL)) != -1 )
    switch (c) {
    case 'g':
      pars->in_geno = optarg;
      break;
    case 'p':
      pars->in_probs = true;
      break;
    case 'l':
      pars->in_logscale = true;
      pars->in_probs = true;
      break;
    case 'n':
      pars->n_ind = atoi(optarg);
      break;
    case 's':
      pars->n_sites = atoi(optarg);
      break;
    case 'a':
      pars->in_pos = optarg;
      pars->in_pos_header = false;
      break;
    case 'A':
      pars->in_pos = optarg;
      pars->in_pos_header = true;
      break;
    case 'd':
      pars->max_kb_dist = atoi(optarg);
      break;
    case 'D':
      pars->max_snp_dist = atoi(optarg);
      break;
    case 'f':
      pars->min_maf = atof(optarg);
      break;
    case 'm':
      pars->ignore_miss_data = true;
      break;
    case 'c':
      pars->call_geno = true;
      break;
    case 'N':
      pars->N_thresh = atof(optarg);
      pars->call_geno = true;
      break;
    case 'C':
      pars->call_thresh = atof(optarg);
      pars->call_geno = true;
      break;
    case 'r':
      pars->rnd_sample = atof(optarg);
      break;
    case 'S':
      pars->seed = atoi(optarg);
      break;
    case 'x':
      pars->extend_out = true;
      break;
    case 'o':
      pars->out = optarg;
      break;
      break;
    case 't':
      pars->n_threads = atoi(optarg);
      break;
    case 'V':
      pars->verbose = atoi(optarg);
      break;
    default:
      exit(-1);
    }


  if(pars->verbose >= 1) {
    fprintf(stderr, "==> Input Arguments:\n");
    fprintf(stderr, "\tgeno: %s\n\tprobs: %s\n\tlog_scale: %s\n\tn_ind: %lu\n\tn_sites: %lu\n\tpos: %s (%s header)\n\tmax_kb_dist (kb): %lu\n\tmax_snp_dist: %lu\n\tmin_maf: %f\n\tignore_miss_data: %s\n\tcall_geno: %s\n\tN_thresh: %f\n\tcall_thresh: %f\n\trnd_sample: %f\n\tseed: %lu\n\textend_out: %s\n\tout: %s\n\tn_threads: %d\n\tverbose: %d\n\tversion: %s (%s @ %s)\n\n",
	    pars->in_geno,
	    pars->in_probs ? "true":"false",
	    pars->in_logscale ? "true":"false",
	    pars->n_ind,
	    pars->n_sites,
	    pars->in_pos,
	    pars->in_pos_header ? "WITH" : "WITHOUT",
	    pars->max_kb_dist,
	    pars->max_snp_dist,
	    pars->min_maf,
	    pars->ignore_miss_data ? "true":"false",
	    pars->call_geno ? "true":"false",
	    pars->N_thresh,
	    pars->call_thresh,
	    pars->rnd_sample,
	    pars->seed,
	    pars->extend_out ? "true":"false",
	    pars->out,
	    pars->n_threads,
	    pars->verbose,
	    version, __DATE__, __TIME__);
  }
  if(pars->verbose > 4)
    fprintf(stderr, "==> Verbose values greater than 4 for debugging purpose only. Expect large amounts of info on screen\n");
  


  /////////////////////
  // Check Arguments //
  /////////////////////
  if(pars->in_geno == NULL)
    error(__FUNCTION__, "genotype input file (--geno) missing!");
  if(pars->n_ind == 0)
    error(__FUNCTION__, "number of individuals (--n_ind) missing!");
  if(pars->n_sites == 0)
    error(__FUNCTION__, "number of sites (--n_sites) missing!");
  if(pars->in_pos == NULL && pars->max_kb_dist > 0)
    error(__FUNCTION__, "position file necessary in order to filter by maximum distance!");
  if(pars->min_maf < 0 || pars->min_maf > 1)
    error(__FUNCTION__, "minimum allele frequency must be in [0,1]!");
  if(pars->call_geno && !pars->in_probs)
    error(__FUNCTION__, "can only call genotypes from likelihoods/probabilities!");
  if(pars->rnd_sample <= 0 || pars->rnd_sample > 1)
    error(__FUNCTION__, "proportion of comparisons to sample must be in ]0,1]!");
  if(pars->n_threads < 1)
    error(__FUNCTION__, "number of threads cannot be less than 1!");
}
