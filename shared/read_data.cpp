
#include "read_data.hpp"

/*
Function to read both called genotypes (1 field per site and indiv), genotype lkls or genotype post probs (3 fields per site and indiv)
  in_geno           : file name to read from
  in_bin            : is input file in binary?
  in_probs          : is input GLs or genot posterior probs?
  in_logscale       : are the probs in log scale? (is updated to true at the end because function always returns genotypes in logscale)
  n_ind             : number of individuals
  n_sites           : number of sites
*/
double*** read_geno(char *in_geno, bool in_bin, bool in_probs, bool *in_logscale, uint64_t n_ind, uint64_t n_sites){
  uint64_t n_fields;
  // Depending on input we will have either 1 or 3 genot
  uint64_t n_geno = (in_probs ? N_GENO : 1);
  double *t, *ptr;
  char *buf = init_ptr(BUFF_LEN, (const char*) '\0');

  // Allocate memory
  double ***geno = init_ptr(n_ind, n_sites, N_GENO, -INF);
  
  // Open GENO file
  gzFile in_geno_fh = open_gzfile(in_geno, in_bin ? "rb" : "r");
  if(in_geno_fh == NULL)
    error(__FUNCTION__, "cannot open GENO file!");

  for(uint64_t s = 0; s < n_sites; s++){
    if(in_bin){
      for(uint64_t i = 0; i < n_ind; i++){
	if( gzread(in_geno_fh, geno[i][s], N_GENO * sizeof(double)) != N_GENO * sizeof(double) )
	  error(__FUNCTION__, "cannot read binary GENO file. Check GENO file and number of sites!");
	if(!*in_logscale)
	  conv_space(geno[i][s], N_GENO, log);
	// Normalize GL
	post_prob(geno[i][s], geno[i][s], NULL, N_GENO);
	// Check if OK
        if(isnan(geno[i][s][0]) ||
           isnan(geno[i][s][1]) ||
           isnan(geno[i][s][2]))
          error(__FUNCTION__, "NaN found! Is the file format correct?");
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
      if(!n_fields || (s == 0 && n_fields < n_ind * n_geno)){
	fprintf(stderr, "> Header found! Skipping line...\n");
        if(s != 0){
          warn(__FUNCTION__, " header found but not on first line. Is this an error?");
          fprintf(stderr, "%s/n", buf);
        }
        s--;
	continue;
      }

      if(n_fields < n_ind * n_geno) {
	fprintf(stderr, "\tline: %s\n\tt[0]: %f\n\tt[1]: %f\n", buf, t[0], t[1]);
	fprintf(stderr, "\tn_line: %lu\n\tfields: %lu\n\tn_ind: %lu\n\tn_geno: %lu\n", s, n_fields, n_ind, n_geno);
	error(__FUNCTION__, "wrong GENO file format. Less fields than expected!");
      }
      
      // Use last "n_ind * n_geno" columns
      ptr = t + (n_fields - n_ind * n_geno);
      
      for(uint64_t i = 0; i < n_ind; i++){
	if(in_probs){
          for(uint64_t g = 0; g < N_GENO; g++)
            geno[i][s][g] = (*in_logscale ? ptr[i*N_GENO+g] : log(ptr[i*N_GENO+g]));
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
  // Read_geno always returns genos in logscale
  *in_logscale = true;
  return geno;
}



int read_split(char *in_file, char*** out, uint64_t* n_rows, uint64_t* n_cols, const char* sep){
  char** buf;

  // Read file
  *n_cols = 0;
  *n_rows = read_file(in_file, &buf);
  if(buf == NULL)
    error(__FUNCTION__, "cannot open file!");

  // Split string into fields
  for(uint64_t i = 0; i < *n_rows; i++) {
    uint64_t n = split(buf[i], sep, &out[i]);
    // If first field, save number of fields
    if(*n_cols == 0)
      *n_cols = n;
    // Check if number if fields is always the same
    if(*n_cols != n)
      error(__FUNCTION__, "invalid number of fields in file!");
  }

  // Clean-up
  delete [] buf;
  return 0;
}



double* read_dist(char *in_pos, uint64_t n_sites){
  uint64_t n_fields = 0;

  // Allocate memory for POS
  double *pos_dist = init_ptr(n_sites, (double) INFINITY);
  // Allocate memory for FILE
  char ***buf = init_ptr(n_sites, 0, 0, (const char*) '\0');

  // Read file
  read_split(in_pos, buf, &n_sites, &n_fields);
  if(buf == NULL)
    error(__FUNCTION__, "cannot open file!");
  if(n_fields < 2)
    error(__FUNCTION__, "wrong POS file format!");

  char *prev_chr = init_ptr(BUFF_LEN, (const char*) '\0');
  uint64_t prev_pos = 0;

  for(uint64_t s = 0; s < n_sites; s++){
    // Check if header and skip
    if(strtod(buf[s][1], NULL) == 0){
      fprintf(stderr, "> Header found! Skipping line...\n");
      if(s != 0){
	warn(__FUNCTION__, " header found but not on first line. Is this an error?");
	fprintf(stderr, "%s\t%s/n", buf[s][0], buf[s][1]);
      }
      s--;
      continue;
    }

    // If first CHR to be parsed, store it
    if(strlen(prev_chr) == 0)
      strcpy(prev_chr, buf[s][0]);
    
    // Check if CHR changed
    if(strcmp(prev_chr, buf[s][0]) == 0){
      pos_dist[s] = strtod(buf[s][1], NULL) - prev_pos;
      if(pos_dist[s] < 1)
	error(__FUNCTION__, "invalid distance between adjacent sites!");
    }else{
      pos_dist[s] = INFINITY;
      strcpy(prev_chr, buf[s][0]);
    }
    prev_pos = strtoul(buf[s][1], NULL, 0);
  }

  delete [] buf;
  delete [] prev_chr;

  return pos_dist;
}
