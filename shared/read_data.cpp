
#include "read_data.hpp"
#include "gen_func.hpp"

// Reads both called genotypes (1 field per site and indiv), genotype lkls or genotype post probs (3 fields per site and indiv)
double*** read_geno(char *in_geno, bool in_bin, bool in_probs, bool in_logscale, uint64_t n_ind, uint64_t n_sites){
  uint64_t n_fields;
  // Depending on input we will have either 1 or 3 genot
  uint64_t n_geno = (in_probs ? N_GENO : 1);
  double *t;
  double *ptr;
  char *buf = init_ptr(BUFF_LEN, (const char*) '\0');

  // Allocate memory
  double ***geno = init_ptr(n_ind, n_sites+1, N_GENO, -INF);
  
  // Open GENO file
  gzFile in_geno_fh = open_gzfile(in_geno, in_bin ? "rb" : "r");
  if(in_geno_fh == NULL)
    error(__FUNCTION__, "cannot open GENO file!");

  for(uint64_t s = 1; s <= n_sites; s++){
    if(in_bin){
      for(uint64_t i = 0; i < n_ind; i++){
	if( gzread(in_geno_fh, geno[i][s], N_GENO * sizeof(double)) != N_GENO * sizeof(double) )
	  error(__FUNCTION__, "cannot read binary GENO file. Check GENO file and number of sites!");
	if(!in_logscale)
	  conv_space(geno[i][s], N_GENO, log);
	// Normalize GL
	post_prob(geno[i][s], geno[i][s], NULL, N_GENO);
      }
    }
    else{
      if( gzgets(in_geno_fh, buf, BUFF_LEN) == NULL)
	error(__FUNCTION__, "cannot read GZip GENO file. Check GENO file and number of sites!");
      // Remove trailing newline
      chomp(buf);
      // Check if empty
      if(strlen(buf) == 0)
	continue;
      // Parse input line into array
      n_fields = split(buf, (const char*) " \t", &t);

      // Check if header and skip
      if(!n_fields){
	s--;
	fprintf(stderr, "> Header found! Skipping line...\n");
	continue;
      }

      if(n_fields < n_ind * n_geno)
	error(__FUNCTION__, "wrong GENO file format. Less fields than expected!");
      
      // Use last "n_ind * n_geno" columns
      ptr = t + (n_fields - n_ind * n_geno);
      
      for(uint64_t i = 0; i < n_ind; i++){
	if(in_probs){
          for(uint64_t g = 0; g < N_GENO; g++)
            geno[i][s][g] = (in_logscale ? ptr[i*N_GENO+g] : log(ptr[i*N_GENO+g]));
	}else{
	  int g = (int) ptr[i];
	  if(g >= 0){
	    if(g > 2)
	      error(__FUNCTION__, "wrong GENO file format. Genotypes must be coded as {-1,0,1,2} !");
	    geno[i][s][g] = log(1);
	  }else
	    geno[i][s][0] = geno[i][s][1] = geno[i][s][2] = log((double) 1/N_GENO);
        }
	
	// Normalize GL
	post_prob(geno[i][s], geno[i][s], NULL, N_GENO);
      }

      // Free memory
      delete [] t;
    }
  }

  // Check if file is at EOF
  gzread(in_geno_fh, buf, 1);
  if(!gzeof(in_geno_fh))
    error(__FUNCTION__, "GENO file not at EOF. Check GENO file and number of sites!");

  gzclose(in_geno_fh);
  delete [] buf;
  return geno;
}







double* read_pos(char *in_pos, uint64_t n_sites){
  uint64_t n_fields;
  char **t;
  char *buf = init_ptr(BUFF_LEN, (const char*) '\0');

  char *prev_chr = init_ptr(BUFF_LEN, (const char*) '\0');
  uint64_t prev_pos = 0;

  // Allocate memory
  double *pos_dist = init_ptr(n_sites+1, (double) INFINITY);

  // Open file
  gzFile in_pos_fh = open_gzfile(in_pos, "r");
  if(in_pos_fh == NULL)
    error(__FUNCTION__, "cannot open POS file!");

  for(uint64_t s = 1; s <= n_sites; s++){
    if( gzgets(in_pos_fh, buf, BUFF_LEN) == NULL)
      error(__FUNCTION__, "cannot read next site from POS file!");
    // Remove trailing newline
    chomp(buf);
    // Check if empty
    if(strlen(buf) == 0)
      continue;
    // Parse input line into array
    n_fields = split(buf, (const char*) " \t", &t);

    // Check if header and skip
    if(!n_fields || strtod(t[1], NULL)==0){
      s--;
      fprintf(stderr, "> Header found! Skipping line...\n");
      continue;
    }

    if(n_fields < 2)
      error(__FUNCTION__, "wrong POS file format!");

    // If first chr to be parsed
    if(strlen(prev_chr) == 0)
      strcpy(prev_chr, t[0]);
    
    if(strcmp(prev_chr, t[0]) == 0){
      pos_dist[s] = strtod(t[1], NULL) - prev_pos;
      if(pos_dist[s] < 1)
	error(__FUNCTION__, "invalid distance between adjacent sites!");
    }else{
      pos_dist[s] = INFINITY;
      strcpy(prev_chr, t[0]);
    }
    prev_pos = strtoul(t[1], NULL, 0);

    // Free memory
    for(uint64_t cnt = 0; cnt < n_fields; cnt++)
      delete [] t[cnt];
    delete [] t;
  }

  // Check if file is at EOF
  gzread(in_pos_fh, buf, 1);
  if(!gzeof(in_pos_fh))
    error(__FUNCTION__, "POS file not at EOF. Check POS file and number of sites!");

  gzclose(in_pos_fh);
  delete [] buf;
  delete [] prev_chr;

  return pos_dist;
}
