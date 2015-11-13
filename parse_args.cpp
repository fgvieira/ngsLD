
#include <getopt.h>
#include "ngsLD.hpp"


void init_pars(params *pars) {
  pars->in_geno = NULL;
  pars->in_bin = false;
  pars->in_probs = false;
  pars->in_logscale = false;
  pars->n_ind = 0;
  pars->n_sites = 0;
  pars->pos = NULL;
  pars->max_dist = 100;
  pars->min_maf = 1e-3;
  pars->call_geno = false;
  pars->N_thresh = 0;
  pars->call_thresh = 0;
  pars->rnd_sample = 1;
  pars->seed = time(NULL) + rand() % 1000;
  pars->out = NULL;
  pars->out_fh = stdout;
  pars->n_threads = 1;
  pars->version = false;
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
      {"pos", required_argument, NULL, 'Z'},
      {"max_dist", required_argument, NULL, 'd'},
      {"min_maf", required_argument, NULL, 'f'},
      {"call_geno", no_argument, NULL, 'c'},
      {"N_thresh", required_argument, NULL, 'N'},
      {"call_thresh", required_argument, NULL, 'C'},
      {"rnd_sample", required_argument, NULL, 'r'},
      {"seed", required_argument, NULL, 'S'},
      {"out", required_argument, NULL, 'o'},
      {"n_threads", required_argument, NULL, 'x'},
      {"version", no_argument, NULL, 'v'},
      {"verbose", required_argument, NULL, 'V'},
      {0, 0, 0, 0}
    };
  
  int c = 0;
  while ( (c = getopt_long_only(argc, argv, "g:pln:s:Z:d:f:cN:C:r:S:o:x:vV:", long_options, NULL)) != -1 )
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
    case 'Z':
      pars->pos = optarg;
      break;
    case 'd':
      pars->max_dist = atof(optarg);
      break;
    case 'f':
      pars->min_maf = atof(optarg);
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
    case 'o':
      pars->out = optarg;
      break;
    case 'x':
      pars->n_threads = atoi(optarg);
      break;
    case 'v':
      pars->version = true;
      break;
    case 'V':
      pars->verbose = atoi(optarg);
      break;
    default:
      exit(-1);
    }


  if(pars->verbose >= 1) {
    fprintf(stderr, "==> Input Arguments:\n");
    fprintf(stderr, "\tgeno: %s\n\tprobs: %s\n\tlog_scale: %s\n\tn_ind: %lu\n\tn_sites: %lu\n\tpos: %s\n\tmax_dist (kb): %lu\n\tmin_maf: %f\n\tcall_geno: %s\n\tN_thresh: %f\n\tcall_thresh: %f\n\trnd_sample: %f\n\tseed: %lu\n\tout: %s\n\tn_threads: %d\n\tversion: %s\n\tverbose: %d\n\n",
	    pars->in_geno,
	    pars->in_probs ? "true":"false",
	    pars->in_logscale ? "true":"false",
	    pars->n_ind,
	    pars->n_sites,
	    pars->pos,
	    pars->max_dist,
	    pars->min_maf,
	    pars->call_geno ? "true":"false",
	    pars->N_thresh,
	    pars->call_thresh,
	    pars->rnd_sample,
	    pars->seed,
	    pars->out,
	    pars->n_threads,
	    pars->version ? "true":"false",
	    pars->verbose
	    );
  }
  if(pars->verbose > 4)
    fprintf(stderr, "==> Verbose values greater than 4 for debugging purpose only. Expect large amounts of info on screen\n");
  


  /////////////////////
  // Check Arguments //
  /////////////////////
  if(pars->version) {
    fprintf(stderr, "ngsLD v%s\nCompiled on %s @ %s\n", version, __DATE__, __TIME__);
    exit(0);
  }
  if(pars->in_geno == NULL)
    error(__FUNCTION__, "genotype input file (--geno) missing!");
  if(pars->n_ind == 0)
    error(__FUNCTION__, "number of individuals (--n_ind) missing!");
  if(pars->n_sites == 0)
    error(__FUNCTION__, "number of sites (--n_sites) missing!");
  if(pars->pos == NULL && pars->max_dist > 0)
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
