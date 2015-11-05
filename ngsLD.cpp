/*
 *
 * ngsLD - NGS data individual inbreeding coefficients estimation.
 * Copyright (C) 2012  Filipe G. Vieira
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include<pthread.h>
#include "ngsLD.hpp"

char const* version = "0.0.1b";


int main (int argc, char** argv) {
  /////////////////////
  // Parse Arguments //
  /////////////////////
  params* pars = new params;
  init_pars(pars);
  parse_cmd_args(pars, argc, argv);

  if(pars->version) {
    fprintf(stderr, "ngsLD v%s\nCompiled on %s @ %s", version, __DATE__, __TIME__);
    exit(0);
  }

  
  
  ///////////////////////
  // Check input files //
  ///////////////////////
  // Get file total size
  struct stat st;
  if(stat(pars->in_geno, &st) != 0)
    error(__FUNCTION__, "cannot check GENO file size!");

  if(strcmp(strrchr(pars->in_geno, '.'), ".gz") == 0){
    if(pars->verbose >= 1)
      fprintf(stderr, "==> GZIP input file (not BINARY)\n");
    pars->in_bin = false;
  }else if(pars->n_sites == st.st_size/sizeof(double)/pars->n_ind/N_GENO){
    if(pars->verbose >= 1)
      fprintf(stderr, "==> BINARY input file (always lkl)\n");
    pars->in_bin = true;
    pars->in_probs = true;
  }else
    error(__FUNCTION__, "invalid/corrupt genotype input file!");



  /////////////////////////////////////////////
  // Declare and initialize output variables //
  /////////////////////////////////////////////
  pars->expected_geno = init_ptr(pars->n_sites+1, pars->n_ind, (double) -1);
  pars->geno_lkl = init_ptr(pars->n_sites+1, pars->n_ind * N_GENO, (double) -1);
  pth_struct **pth = new pth_struct*[pars->n_sites+1];
  // Open filehandle
  if(pars->out != NULL)
    pars->out_fh = fopen(pars->out, "w");
  if(pars->out_fh == NULL)
    error(__FUNCTION__, "cannot open output file!");
  //fprintf(pars->out_fh, "#site1\tsite2\tdist\tr^2_ExpG\tD\tD'\tr^2\n");



  /////////////////////
  // Read input data //
  /////////////////////
  // Read data from GENO file
  if(pars->verbose >= 1)
    fprintf(stderr, "> Reading data from file...\n");
  pars->in_geno_lkl = read_geno(pars->in_geno, pars->in_bin, pars->in_probs, pars->in_logscale, pars->n_ind, pars->n_sites);

  // Read positions from file
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Getting sites coordinates\n");
  if(pars->pos)
    pars->pos_dist = read_pos(pars->pos, pars->n_sites);
  else
    pars->pos_dist = init_ptr(pars->n_sites+1, INFINITY);

  // Read labels
  if(read_file(pars->pos, &pars->labels) != pars->n_sites)
    error(__FUNCTION__, "number of labels does not match number of sites!");
  // Fix labels...
  char* ptr;
  for(uint64_t s = 1; s <= pars->n_sites; s++){
    ptr = strchr(pars->labels[s-1], '\t');
    *ptr = ':';
  }

  // Data pre-processing...
  for(uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t s = 1; s <= pars->n_sites; s++){
      // Call genotypes
      if(pars->call_geno)
	call_geno(pars->in_geno_lkl[i][s], N_GENO);

      // Convert to normal space
      conv_space(pars->in_geno_lkl[i][s], N_GENO, exp);

      // Caclulate expected genotypes
      pars->expected_geno[s][i] = pars->in_geno_lkl[i][s][1] + 2*pars->in_geno_lkl[i][s][2];

      // Convert matrix format
      memcpy(&pars->geno_lkl[s][i*N_GENO], pars->in_geno_lkl[i][s], N_GENO * sizeof(double));
    }
  


  //////////////////
  // Analyze Data //
  //////////////////
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Launching threads...\n");

  // Create thread pool
  if( (pars->thread_pool = threadpool_create(pars->n_threads, pars->n_sites, 0)) == NULL )
    error(__FUNCTION__, "failed to create thread pool!");

  for(uint64_t s1 = 1; s1 < pars->n_sites; s1++){
    // Fill in pthread structure
    pth[s1] = new pth_struct;
    pth[s1]->pars = pars;
    pth[s1]->site = s1;

    // Add task to thread pool
    int ret = threadpool_add(pars->thread_pool, calc_pair_LD, (void*) pth[s1], 0);
    if(ret == -1)
      error(__FUNCTION__, "invalid thread pool!");
    else if(ret == -2)
      error(__FUNCTION__, "thread pool lock failure!");
    else if(ret == -3)
      error(__FUNCTION__, "queue full!");
    else if(ret == -4)
      error(__FUNCTION__, "thread pool is shutting down!");
    else if(ret == -5)
      error(__FUNCTION__, "thread failure!");
  }



  //////////////////////////
  // Wait for all threads //
  //////////////////////////
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Waiting for all threads to finish...\n");
  
  threadpool_wait(pars->thread_pool, 0.1);
  if(threadpool_destroy(pars->thread_pool, threadpool_graceful) != 0)
    error(__FUNCTION__, "cannot free thread pool!");



  /////////////////
  // Free Memory //
  /////////////////
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Freeing memory...\n");

  fclose(pars->out_fh);
  //free_ptr((void*) pars->in_geno);
  free_ptr((void***) pars->in_geno_lkl, pars->n_ind, pars->n_sites+1);
  free_ptr((void**) pars->geno_lkl, pars->n_sites+1);
  free_ptr((void**) pars->expected_geno, pars->n_sites+1);
  free_ptr((void**) pars->labels, pars->n_sites);
  free_ptr((void*) pars->pos_dist);
  free_ptr((void**) pth);

  if(pars->verbose >= 1)
    fprintf(stderr, "Done!\n");
  delete pars;

  return 0;
}





void calc_pair_LD (void *pth){
  pth_struct* p = (pth_struct*) pth;
  uint64_t s1 = p->site;
  uint64_t s2 = s1;
  double dist = 0;
  double LD_GL[3];

  // Calc LD for pairs of SNPs < max_dist
  do{
    s2++;
    dist += p->pars->pos_dist[s2];

    if(p->pars->verbose > 8)
      fprintf(stderr, "%lu\t%s\t%lu\t%s\t%.0f\t%s\t%s\n", s1, p->pars->labels[s1-1], s2, p->pars->labels[s2-1], dist, join(p->pars->expected_geno[s1],p->pars->n_ind,","), join(p->pars->expected_geno[s2],p->pars->n_ind,","));
    
    if(p->pars->max_dist > 0 && dist >= p->pars->max_dist*1000)
      break;

    // Calculate LD using bcftools algorithm
    bcf_pair_LD(LD_GL, p->pars->geno_lkl[s1], p->pars->geno_lkl[s2], p->pars->n_ind);

    /*
    if( isnan(LD_GL[1]) && isnan(LD_GL[2]) ){
      //fprintf(p->pars->out_fh, "%s\t%s\t%.0f\n", p->pars->labels[s1-1], p->pars->labels[s2-1], dist);
      continue;
    }
    */

    fprintf(p->pars->out_fh, "%s\t%s\t%.0f\t%f\t%f\t%f\t%f\n",
	    p->pars->labels[s1-1],
	    p->pars->labels[s2-1],
	    dist,
	    pearson_r(p->pars->expected_geno[s1], p->pars->expected_geno[s2], p->pars->n_ind),
	    LD_GL[0],
	    LD_GL[1],
	    LD_GL[2]
	    );
  } while (s2 < p->pars->n_sites);

  delete p;
}





double pearson_r (double *s1, double *s2, uint64_t n_ind){
  return pow(gsl_stats_correlation(s1, 1, s2, 1, n_ind), 2);
}



#define ITER_MAX 50
#define EPS 1e-5


// Adapted from BCFTOOLS:
// https://github.com/lh3/samtools/blob/6bbe1609e10f27796e5bf29ac3207bb2e35ceac8/bcftools/em.c#L266-L310
void bcf_pair_LD (double LD[3], double *s1, double *s2, uint64_t n_ind)
{
  double f[4], flast[4], f0[2];

  // set the initial value
  f[0] = f[1] = f[2] = f[3] = -1.0;
  f0[0] = est_freq(n_ind, s1);
  f0[1] = est_freq(n_ind, s2);
  f[0] = (1 - f0[0]) * (1 - f0[1]); f[3] = f0[0] * f0[1];
  f[1] = (1 - f0[0]) * f0[1]; f[2] = f0[0] * (1 - f0[1]);

  // iteration
  for(uint64_t i = 0; i < ITER_MAX; i++) {
    double eps = 0;
    memcpy(flast, f, 4 * sizeof(double));
    pair_freq_iter(n_ind, s1, s2, f);
    for (uint64_t j = 0; j < 4; j++) {
      double x = fabs(f[j] - flast[j]);
      if (x > eps) eps = x;
    }

    if(eps < EPS)
      break;
  }

  // calculate r^2
  double r2, D, Dp, p[2];
  p[0] = f[0] + f[1];
  p[1] = f[0] + f[2];
  D = f[0] * f[3] - f[1] * f[2]; // P_AB * P_BA - P_AA * P_BB
  // or
  //D = f[0] - p[0] * p[1];      // P_AB - P_A * P_B
  Dp = D / (D < 0 ? min(p[0]*p[1], (1-p[0])*(1-p[1])) : min(p[0]*(1-p[1]), (1-p[0])*p[1]) );
  r2 = pow(D / sqrt(p[0] * p[1] * (1-p[0]) * (1-p[1])), 2);

  /*
  if(isnan(r2))
    r2 = -1.0;
  */

  LD[0] = D;
  LD[1] = Dp;
  LD[2] = r2;
}



#define _G1(h, k) ((h>>1&1) + (k>>1&1))
#define _G2(h, k) ((h&1) + (k&1))

// 0: the previous site; 1: the current site
int pair_freq_iter(int n, double *s1, double *s2, double f[4])
{
  double ff[4];
  int i, k, h;
  //fprintf(stderr, "%lf,%lf,%lf,%lf\n", f[0], f[1], f[2], f[3]);
  memset(ff, 0, 4 * sizeof(double));
  for (i = 0; i < n; ++i) {
    double *p[2], sum, tmp;
    p[0] = s1 + i * 3;
    p[1] = s2 + i * 3;
    for (k = 0, sum = 0.; k < 4; ++k)
      for (h = 0; h < 4; ++h)
	sum += f[k] * f[h] * p[0][_G1(k,h)] * p[1][_G2(k,h)];
    for (k = 0; k < 4; ++k) {
      tmp = f[0] * (p[0][_G1(0,k)] * p[1][_G2(0,k)] + p[0][_G1(k,0)] * p[1][_G2(k,0)])
	+ f[1] * (p[0][_G1(1,k)] * p[1][_G2(1,k)] + p[0][_G1(k,1)] * p[1][_G2(k,1)])
	+ f[2] * (p[0][_G1(2,k)] * p[1][_G2(2,k)] + p[0][_G1(k,2)] * p[1][_G2(k,2)])
	+ f[3] * (p[0][_G1(3,k)] * p[1][_G2(3,k)] + p[0][_G1(k,3)] * p[1][_G2(k,3)]);
      ff[k] += f[k] * tmp / sum;
    }
  }
  for (k = 0; k < 4; ++k) f[k] = ff[k] / (2 * n);
  return 0;
}



// estimate site allele frequency in a very naive and inaccurate way
double est_freq(int n, const double *pdg)
{
  int i, gcnt[3], tmp1;
  // get a rough estimate of the genotype frequency
  gcnt[0] = gcnt[1] = gcnt[2] = 0;
  for (i = 0; i < n; ++i) {
    const double *p = pdg + i * 3;
    if (p[0] != 1. || p[1] != 1. || p[2] != 1.) {
      int which = p[0] > p[1]? 0 : 1;
      which = p[which] > p[2]? which : 2;
      ++gcnt[which];
    }
  }
  tmp1 = gcnt[0] + gcnt[1] + gcnt[2];
  return (tmp1 == 0)? -1.0 : (.5 * gcnt[1] + gcnt[2]) / tmp1;
}
