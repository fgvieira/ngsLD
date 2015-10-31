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
 
    if(pars->verbose > 5)
      printf("Launching thread for site %lu!\n", s1);

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
    error(__FUNCTION__, "number of labels does not match number of individuals!");

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

    if(p->pars->verbose > 7)
      printf("\t%lu <=> %lu: %f\n", s1, s2, dist);

    p->pars->matrixLD[s1][s2] = pearson_r(p->pars->expected_geno[s1], p->pars->expected_geno[s2], p->pars->n_ind);
    p->pars->matrixLD[s2][s1] = dist;
  } while (s2 < p->pars->n_sites);

  delete p;
}





double pearson_r (double *s1, double *s2, uint64_t n_ind){
  return gsl_stats_correlation(s1, 1, s2, 1, n_ind);
}


/*
double bcftools (){
 https://github.com/lh3/samtools/blob/6bbe1609e10f27796e5bf29ac3207bb2e35ceac8/bcftools/em.c
}
*/
