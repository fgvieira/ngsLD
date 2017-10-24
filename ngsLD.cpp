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
  }else{
    if(pars->verbose >= 1)
      fprintf(stderr, "==> BINARY input file (always lkl)\n");
    pars->in_bin = true;
    pars->in_probs = true;

    if(pars->n_sites != st.st_size/sizeof(double)/pars->n_ind/N_GENO)
      error(__FUNCTION__, "invalid/corrupt genotype input file!");
  }



  ////////////////////
  // Prepare output //
  ////////////////////
  // Initialize MAF array
  pars->maf = init_ptr(pars->n_sites+1, (double) -1);
  // Initialize output variables
  pars->expected_geno = init_ptr(pars->n_sites+1, pars->n_ind, (double) -1);
  // Initialize random seed generator
  gsl_rng* rnd_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rnd_gen, pars->seed);
  
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
  double ***tmp = read_geno(pars->in_geno, pars->in_bin, pars->in_probs, &pars->in_logscale, pars->n_ind, pars->n_sites);
  pars->geno_lkl = transp_matrix(tmp, pars->n_ind, pars->n_sites+1);
  free_ptr((void**) tmp, pars->n_ind);

  // Call genotypes
  if(pars->call_geno){
    if(pars->verbose >= 1)
      fprintf(stderr, "> Calling genotypes...\n");
    for(uint64_t i = 0; i < pars->n_ind; i++)
      for(uint64_t s = 1; s <= pars->n_sites; s++)
	call_geno(pars->geno_lkl[s][i], N_GENO, pars->in_logscale, pars->N_thresh, pars->call_thresh, 0);
  }

  // Calculate MAF
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Calculating MAF for all sites...\n");
  for(uint64_t s = 1; s <= pars->n_sites; s++)
    pars->maf[s] = est_maf(pars->n_ind, pars->geno_lkl[s], (double*) NULL);
  
  // Data pre-processing...
  for(uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t s = 1; s <= pars->n_sites; s++){
      // Convert to normal space (since GENOs are in LOG from read_geno)
      conv_space(pars->geno_lkl[s][i], N_GENO, exp);

      // Calculate expected genotypes
      pars->expected_geno[s][i] = pars->geno_lkl[s][i][1] + 2*pars->geno_lkl[s][i][2];
    }


  // Read positions from file
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Getting sites coordinates\n");
  if(pars->pos)
    pars->pos_dist = read_pos(pars->pos, pars->n_sites);
  else
    pars->pos_dist = init_ptr(pars->n_sites+1, (double) INFINITY);

  // Read labels
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Getting sites labels\n");
  pars->labels = read_file(pars->pos, 0, pars->n_sites, BUFF_LEN);
  if(pars->labels == NULL)
    error(__FUNCTION__, "number of labels does not match number of sites!");
  // Fix labels...
  char* ptr;
  for(uint64_t s = 1; s <= pars->n_sites; s++){
    ptr = strchr(pars->labels[s-1], '\t');
    if(ptr != NULL)
      *ptr = ':';
  }



  // DEBUG
  if(pars->verbose >= 7)
    for(uint64_t s = 1; s <= min(10, pars->n_sites); s++)
      fprintf(stderr, "%lu\t%s\t%f (%f %f %f)\n", s, pars->labels[s-1], pars->maf[s], pars->geno_lkl[s][0][0], pars->geno_lkl[s][0][1], pars->geno_lkl[s][0][2]);



  //////////////////
  // Analyze Data //
  //////////////////
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Launching threads...\n");

  // Create thread pool
  if( (pars->thread_pool = threadpool_create(pars->n_threads, pars->n_sites, 0)) == NULL )
    error(__FUNCTION__, "failed to create thread pool!");
  // Allocate memory for array of pthread structures
  pth_struct **pth = new pth_struct*[pars->n_sites+1];

  for(uint64_t s1 = 1; s1 <= pars->n_sites; s1++){
    // Fill in pthread structure
    pth[s1] = new pth_struct;
    pth[s1]->pars = pars;
    pth[s1]->site = s1;
    // Since ngsLD is threaded, in order for the results to be replicable, each thread has its own random number generator.
    pth[s1]->rnd_gen = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(pth[s1]->rnd_gen, draw_rnd(rnd_gen, 0, INF));

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

    if(pars->verbose >= 3)
      if(s1 == 1 ||
         s1 == pars->n_sites ||
         s1 % 100000 == 0)
        fprintf(stderr, "> Launched thread for site %lu (%s)\n", s1, pars->labels[s1-1]);
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
  free_ptr((void***) pars->geno_lkl, pars->n_sites+1, pars->n_ind);
  free_ptr((void**) pars->expected_geno, pars->n_sites+1);
  free_ptr((void**) pars->labels, pars->n_sites);
  free_ptr((void*) pars->pos_dist);
  free_ptr((void*) pars->maf);
  free_ptr((void**) pth); // Only the **pth is freed since the others are freed by each thread.
  gsl_rng_free(rnd_gen);

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
  double LD_GL[4];


  // Calc LD for pairs of SNPs < max_kb_dist
  while (s2 < p->pars->n_sites){
    s2++;
    dist += p->pars->pos_dist[s2];

    if(p->pars->verbose > 8)
      fprintf(stderr, "%lu\t%s\t%lu\t%s\t%.0f\t%f\t%f\t%s\t%s\n", s1, p->pars->labels[s1-1], s2, p->pars->labels[s2-1], dist, p->pars->maf[s1], p->pars->maf[s2], join(p->pars->expected_geno[s1],p->pars->n_ind,","), join(p->pars->expected_geno[s2],p->pars->n_ind,","));

    // Stop if current distance is greater than max_kb_dist
    if(p->pars->max_kb_dist > 0 && p->pars->max_kb_dist*1000 < dist)
      break;
    // Stop if current SNP is greater than max_snp_dist
    if(p->pars->max_snp_dist > 0 && p->pars->max_snp_dist < s2 - s1)
      break;

    // Random sampling
    if(draw_rnd(p->rnd_gen, 0, 1) > p->pars->rnd_sample)
      continue;

    // If either site 1 or site 2 is monomorphic, set LD to "NAN"
    if(p->pars->maf[s1] < p->pars->min_maf || p->pars->maf[s2] < p->pars->min_maf){
      LD_GL[0] = LD_GL[1] = LD_GL[2] = LD_GL[3] = NAN;
    }else{
      // Calculate LD using bcftools algorithm
      LD_GL[0] = pearson_r(p->pars->expected_geno[s1], p->pars->expected_geno[s2], p->pars->n_ind);

      // Calculate LD using bcftools algorithm
      bcf_pair_LD(LD_GL+1, p->pars->geno_lkl[s1], p->pars->geno_lkl[s2], p->pars->maf[s1], p->pars->maf[s2], p->pars->n_ind, false);
    }

    // If sampling, also skip nan
    if(p->pars->rnd_sample != 1 && isnan(LD_GL[0]))
      continue;

    fprintf(p->pars->out_fh, "%s\t%s\t%.0f\t%f\t%f\t%f\t%f\n",
	    p->pars->labels[s1-1],
	    p->pars->labels[s2-1],
	    dist,
	    LD_GL[0],
	    LD_GL[1],
	    LD_GL[2],
	    LD_GL[3]
	    );
  }

  // Free memory
  gsl_rng_free(p->rnd_gen);
  delete p;
}





double pearson_r (double *s1, double *s2, uint64_t n_ind){
  return pow(gsl_stats_correlation(s1, 1, s2, 1, n_ind), 2);
}
