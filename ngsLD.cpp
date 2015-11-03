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
    printf("ngsLD v%s\nCompiled on %s @ %s", version, __DATE__, __TIME__);
    exit(0);
  }

  
  
  ///////////////////////
  // Check input files //
  ///////////////////////
  // Get file total size
  struct stat st;
  if(stat(pars->in_geno, &st) != 0)
    error(__FUNCTION__, "cannot check file size!");

  if(strcmp(strrchr(pars->in_geno, '.'), ".gz") == 0){
    if(pars->verbose >= 1)
      printf("==> GZIP input file (not BINARY)\n");
    pars->in_bin = false;
  }else if(pars->n_sites == st.st_size/sizeof(double)/pars->n_ind/N_GENO){
    if(pars->verbose >= 1)
      printf("==> BINARY input file (always lkl)\n");
    pars->in_bin = true;
    pars->in_probs = true;
  }else
    error(__FUNCTION__, "invalid/corrupt genotype input file!");



  /////////////////////////////////////////////
  // Declare and initialize output variables //
  /////////////////////////////////////////////
  pars->expected_geno = init_ptr(pars->n_sites+1, pars->n_ind, (double) -1);
  pars->matrixLD = init_ptr(pars->n_sites+1, pars->n_sites+1, (double) 0);
  pth_struct **pth = new pth_struct*[pars->n_sites+1];



  /////////////////////
  // Read input data //
  /////////////////////
  // Read data from GENO file
  if(pars->verbose >= 1)
    printf("> Reading data from file...\n");
  pars->geno_lkl = read_geno(pars->in_geno, pars->in_bin, pars->in_probs, pars->in_logscale, pars->n_ind, pars->n_sites);

  // Read positions from file
  if(pars->verbose >= 1)
    printf("==> Getting sites coordinates\n");
  if(pars->pos)
    pars->pos_dist = read_pos(pars->pos, pars->n_sites);
  else
    pars->pos_dist = init_ptr(pars->n_sites+1, INFINITY);
  // Convert position distances to Kb
  for(uint64_t s = 1; s <= pars->n_sites; s++)
    pars->pos_dist[s] /= 1e3;

  // Data pre-processing...
  for(uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t s = 1; s <= pars->n_sites; s++){
      // Call genotypes
      if(pars->call_geno)
	call_geno(pars->geno_lkl[i][s], N_GENO);

      // Convert to normal space
      conv_space(pars->geno_lkl[i][s], N_GENO, exp);

      // Caclulate expected genotypes
      pars->expected_geno[s][i] = pars->geno_lkl[i][s][1] + 2*pars->geno_lkl[i][s][2];
    }
  


  //////////////////
  // Analyze Data //
  //////////////////
  if(pars->verbose >= 1)
    printf("==> Launching threads...\n");

  // Create thread pool
  if( (pars->thread_pool = threadpool_create(pars->n_threads, pars->n_sites, 0)) == NULL )
    error(__FUNCTION__, "failed to create thread pool!");

  for(uint64_t s1 = 1; s1 < pars->n_sites; s1++){
    // Fill in pthread structure
    pth[s1] = new pth_struct;
    pth[s1]->pars = pars;
    pth[s1]->site = s1;
 
    //calc_pair_LD((void*) pth[s1]);

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
    printf("==> Waiting for all threads to finish...\n");
  
  threadpool_wait(pars->thread_pool, 0.1);
  if(threadpool_destroy(pars->thread_pool, threadpool_graceful) != 0)
    error(__FUNCTION__, "cannot free thread pool!");



  //////////////////
  // Print Matrix //
  //////////////////
  if(pars->verbose >= 1)
    printf("==> Printing LD results...\n");

  char **labels;
  if(read_file(pars->pos, &labels) != pars->n_sites)
    error(__FUNCTION__, "number of labels does not match number of sites!");

  FILE* out_fh = fopen(pars->out, "w");
  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open output file!");

  for(uint64_t s1 = 1; s1 <= pars->n_sites; s1++)
    for(uint64_t s2 = s1+1; s2 <= pars->n_sites; s2++)
      if(pars->matrixLD[s1][s2] >= pars->min_r2 && pars->matrixLD[s2][s1] > 0)
	fprintf(out_fh, "%s\t%s\t%f\t%f\n", labels[s1-1], labels[s2-1], pars->matrixLD[s1][s2], pars->matrixLD[s2][s1]);

  fclose(out_fh);



  /////////////////
  // Free Memory //
  /////////////////
  if(pars->verbose >= 1)
    printf("==> Freeing memory...\n");

  // pars struct
  //free_ptr((void*) pars->in_geno);
  free_ptr((void***) pars->geno_lkl, pars->n_ind, pars->n_sites+1);
  free_ptr((void*) pars->pos_dist);
  free_ptr((void**) pth);
  free_ptr((void**) pars->expected_geno, pars->n_sites+1);
  free_ptr((void**) pars->matrixLD, pars->n_sites+1);
  free_ptr((void**) labels, pars->n_sites);
  

  if(pars->verbose >= 1)
    printf("Done!\n");
  delete pars;

  return 0;
}





void calc_pair_LD (void *pth){
  pth_struct* p = (pth_struct*) pth;
  uint64_t s1 = p->site;
  uint64_t s2 = s1;
  double dist = 0;

  // Calc LD for pairs of SNPs < max_dist
  do{
    s2++;
    dist += p->pars->pos_dist[s2];
    
    if(p->pars->max_dist != -1 && dist >= p->pars->max_dist)
      break;

    if(1){
      p->pars->matrixLD[s1][s2] = pearson_r(p->pars->expected_geno[s1], p->pars->expected_geno[s2], p->pars->n_ind);
      //    }else{
      //      p->pars->matrixLD[s1][s2] = bcf_pair_LD(p->pars->geno_lkl[s1], p->pars->geno_lkl[s2], p->pars->n_ind);
    }
    p->pars->matrixLD[s2][s1] = dist;

    if(p->pars->verbose > 7)
      printf("\t%lu <=> %lu: %f (%f)\n", s1, s2, p->pars->matrixLD[s1][s2], dist);
  } while (s2 < p->pars->n_sites);

  delete p;
}





double pearson_r (double *s1, double *s2, uint64_t n_ind){
  return gsl_stats_correlation(s1, 1, s2, 1, n_ind);
}



#define ITER_MAX 50
#define EPS 1e-5


// Adapted from BCFTOOLS:
// https://github.com/lh3/samtools/blob/6bbe1609e10f27796e5bf29ac3207bb2e35ceac8/bcftools/em.c#L266-L310
double bcf_pair_LD (double *s1, double *s2, uint64_t n_ind)
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
  double r, D, p[2];
  p[0] = f[0] + f[1];
  p[1] = f[0] + f[2];
  D = f[0] * f[3] - f[1] * f[2];
  r = sqrt(D * D / (p[0] * p[1] * (1-p[0]) * (1-p[1])));
  //printf("R(%lf,%lf,%lf,%lf)=%lf\n", f[0], f[1], f[2], f[3], r);

  if(isnan(r))
    r = -1.0;
  
  return r;
}



#define _G1(h, k) ((h>>1&1) + (k>>1&1))
#define _G2(h, k) ((h&1) + (k&1))

// 0: the previous site; 1: the current site
int pair_freq_iter(int n, double *s1, double *s2, double f[4])
{
  double ff[4];
  int i, k, h;
  //printf("%lf,%lf,%lf,%lf\n", f[0], f[1], f[2], f[3]);
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
