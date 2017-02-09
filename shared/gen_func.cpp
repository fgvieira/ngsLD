#include "gen_func.hpp"



void warn(const char *func, const char *msg) {
  fflush(stdout);
  fprintf(stderr, "\n=======\nWARNING: [%s] %s\n=======\n\n", func, msg);
  fflush(stderr);
}


void error(const char *func, const char *msg) {
  fflush(stdout);
  fprintf(stderr, "\n=====\nERROR: [%s] %s\n=====\n\n", func, msg);
  perror("\t");
  fflush(stderr);
  exit(-1);
}


bool SIG_COND;
int really_kill = 3;
void handler(int s) {
  if(SIG_COND)
    fprintf(stderr,"\n\"%s\" signal caught! Will try to exit nicely (no more threads are created but we will wait for the current threads to finish).\n", strsignal(s));

  really_kill--;
  if(really_kill > 0)
    fprintf(stderr,"\t-> If you really want to force an unclean exit Ctr+C %d more times\n", really_kill);
  fflush(stderr);

  if(!really_kill)
    exit(0);

  SIG_COND = false;
}


//we are threading so we want make a nice signal handler for ctrl+c
void catch_SIG(){
  struct sigaction sa;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;

  //sigaction(SIGKILL, &sa, 0);
  sigaction(SIGTERM, &sa, 0);
  sigaction(SIGQUIT, &sa, 0);
  //sigaction(SIGABRT, &sa, 0);
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);
}


double check_interv(double value, bool verbose) {
  if (isnan(value))
    error(__FUNCTION__, "value is NaN!\n");

  if(value < EPSILON) {
    value = 0;
    if(verbose && value < 0)
      fprintf(stderr, "\nWARN: value %f < 0!\n", value);
  } else if(value > 1 - EPSILON) {
    value = 1;
    if(verbose && value > 1)
      fprintf(stderr, "\nWARN: value %f > 1!\n", value);
  }

  return value;
}


int array_max_pos(double *array, int size) {
  int res = 0;
  double max = -INFINITY;

  for (int cnt = 0; cnt < size; cnt++) {
    if (array[cnt] > max) {
      res = cnt;
      max = array[cnt];
    }
  }
  return res;
}


int array_min_pos(double *array, int size) {
  int res = 0;
  double min = +INFINITY;

  for (int cnt = 0; cnt < size; cnt++) {
    if (array[cnt] < min) {
      res = cnt;
      min = array[cnt];
    }
  }
  return res;
}


// Transpose matrix
double*** transp_matrix(double ***src, uint64_t A, uint64_t B){
  if(A < 1 || B < 1)
    error(__FUNCTION__, "invalid size of array!");

  double*** dest = init_ptr(B, A, 0, (double)0);
  
  for(uint64_t a = 0; a < A; a++)
    for(uint64_t b = 0; b < B; b++)
      dest[b][a] = src[a][b];

  return dest;
}



double draw_rnd(gsl_rng *r, uint64_t min, uint64_t max) {
  return( min + gsl_rng_uniform(r) * (max - min) );
}



void conv_space(double *geno, int n_geno, double (*func)(double)) {
  for (int g = 0; g < n_geno; g++) {
    geno[g] = func(geno[g]);

    if(geno[g] == -INFINITY)
      geno[g] = -INF;
  }
}



//function does: log(exp(a)+exp(b)) while protecting for underflow
double logsum(double *a, uint64_t n){
  double sum = 0;

  // Find maximum value
  double M = a[0];
  for(uint64_t i = 1; i < n; i++)
    M = max(a[i], M);

  // If all elements -Inf, return -Inf
  if(M == -INFINITY)
    return(-INFINITY);

  for(uint64_t i = 0; i < n; i++)
    sum += exp(a[i] - M);

  return(log(sum) + M);
}


// Special logsum case for size == 2
double logsum(double a, double b){
  double buf[2];
  buf[0] = a;
  buf[1] = b;
  return logsum(buf, 2);
}


// Special logsum case for size == 3
double logsum(double a, double b, double c){
  double buf[3];
  buf[0] = a;
  buf[1] = b;
  buf[2] = c;
  return logsum(buf, 3);
}


// Special logsum case for size == 4
double logsum(double a, double b, double c, double d){
  double buf[4];
  buf[0] = a;
  buf[1] = b;
  buf[2] = c;
  buf[3] = d;
  return logsum(buf, 4);
}



// Remove trailing newlines from strings
void chomp(char *str){
  char *last_char = &str[strlen(str)-1];

  if(strcmp(last_char,"\r") == 0 ||
     strcmp(last_char,"\n") == 0 ||
     strcmp(last_char,"\r\n") == 0 )
    *last_char = '\0';
}



// Open gzFile
gzFile open_gzfile(const char* name, const char* mode, uint64_t buf_size){
  gzFile fh = NULL;

  if( strcmp(name, "-") == 0 )
    fh = gzdopen(fileno(stdin), mode);
  else
    fh = gzopen(name, mode);

  if(fh == NULL)
    return NULL;

  if(gzbuffer(fh, buf_size) < 0)
    return NULL;

  return fh;
}



// Read data from file and place into array
char **read_file(const char *in_file, uint64_t start_line, uint64_t n_lines, uint64_t buff_size){
  uint64_t cnt = 0;
  char buf[buff_size];
  char **ptr = init_ptr(n_lines-start_line, 0, (const char*) '\0');

  // Open file
  gzFile in_file_fh = gzopen(in_file, "r");
  if(in_file_fh == NULL)
    return NULL;

  for(cnt = 0; cnt < n_lines && !gzeof(in_file_fh); cnt++){
    buf[0] = '\0';
    // Read line from file
    gzgets(in_file_fh, buf, buff_size);
    // Remove trailing newline
    chomp(buf);
    // Check if empty
    if(strlen(buf) == 0){
      cnt--;
      continue;
    }
    // Skip until star_line
    if(cnt < start_line){
      cnt--;
      continue;
    }
    // Alloc memory and copy line
    ptr[cnt] = init_ptr(strlen(buf)+1, buf);
  }

  if(cnt != n_lines)
    error(__FUNCTION__, "could not read specified number of lines!");

  gzclose(in_file_fh);
  return ptr;
}



// Read data from file and place into array double
double **read_file(const char *in_file, uint64_t start_line, uint64_t n_lines, int n_cols, uint64_t buff_size){
  char **tmp = read_file(in_file, start_line, n_lines, buff_size);
  if(tmp == NULL)
    return NULL;

  double **ptr = init_ptr(n_lines-start_line, 0, 0.0);
  for(uint64_t cnt = 0; cnt < n_lines; cnt++)
    if(split(tmp[cnt], (const char*) " \t", &ptr[cnt]) != n_cols)
      error(__FUNCTION__, "number of columns do not match!");

  free_ptr((void**) tmp, n_lines);
  return ptr;
}



// New strtok function to allow for empty ("") separators
char *_strtok(char **str, const char *sep){
  size_t pos = 1;
  char *tmp = strdcat(*str, "");

  if(strcmp(sep, "") != 0)
    pos = strcspn(*str, sep);

  *str += pos;
  tmp[pos] = '\0';

  if(strlen(*str) == 0)
    *str = NULL;
  else if(strcmp(sep, "") != 0)
    *str += 1;

  return tmp;
}


/***************************
split()
 function to, given a 
 string (str), splits into 
 array (out) on char (sep)
***************************/
int split(char *str, const char *sep, int **out){
  int i = strlen(str);
  int *buf = init_ptr(i, 0);

  i = 0;
  char *pch;
  char *end_ptr;
  while(str != NULL){
    pch = _strtok(&str, sep);
    if(strlen(pch) == 0){
      delete [] pch;
      continue;
    }
    
    buf[i++] = strtol(pch, &end_ptr, 0);
    // Check if an int
    if(*end_ptr)
      i--;

    delete [] pch;
  }

  *out = init_ptr(i, 0); // FGV: why the need for *out?!?!!?
  memcpy(*out, buf, i*sizeof(int));

  delete [] buf;
  return(i);
}


int split(char *str, const char *sep, float **out){
  int i = strlen(str);
  float *buf = init_ptr(i, (float) 0);

  i = 0;
  char *pch;
  char *end_ptr;
  while(str != NULL){
    pch = _strtok(&str, sep);
    if(strlen(pch) == 0){
      delete [] pch;
      continue;
    }

    buf[i++] = strtof(pch, &end_ptr);
    // Check if float
    if(*end_ptr)
      i--;

    delete [] pch;
  }

  *out = init_ptr(i, (float) 0);
  memcpy(*out, buf, i*sizeof(float));
  delete [] buf;

  return(i);
}


int split(char *str, const char *sep, double **out){
  int i = strlen(str);
  double *buf = init_ptr(i, (double) 0);

  i = 0;
  char *pch;
  char *end_ptr;
  while(str != NULL){
    pch = _strtok(&str, sep);
    if(strlen(pch) == 0){
      delete [] pch;
      continue;
    }

    buf[i++] = strtod(pch, &end_ptr);
    // Check if double
    if(*end_ptr)
      i--;

    delete [] pch;
  }

  *out = init_ptr(i, (double) 0);
  memcpy(*out, buf, i*sizeof(double));
  delete [] buf;

  return(i);
}


int split(char *str, const char *sep, char ***out){
  int i = strlen(str);
  char **buf = init_ptr(i, 0, (char*) NULL);

  i = 0;
  while(str != NULL)
    buf[i++] = _strtok(&str, sep);

  *out = init_ptr(i, 0, (char*) NULL);
  memcpy(*out, buf, i*sizeof(char*));
  delete [] buf;

  return(i);
}



char *join(unsigned short int *array, uint64_t size, const char *sep){
  char *buf = init_ptr(size*10, (char*) NULL);
  uint64_t len = 0;

  sprintf(buf, "%u", array[0]);
  len = strlen(buf);

  for(uint64_t cnt = 1; cnt < size; cnt++){
    sprintf(buf+len, "%s%u", sep, array[cnt]);
    len = strlen(buf);
  }
  
  char *str = init_ptr(len+1, (char*) NULL);
  strcpy(str, buf);
  delete [] buf;

  return str;
}



char *join(uint64_t *array, uint64_t size, const char *sep){
  char *buf = init_ptr(size*25, (char*) NULL);
  uint64_t len = 0;

  sprintf(buf, "%lu", array[0]);
  len = strlen(buf);

  for(uint64_t cnt = 1; cnt < size; cnt++){
    sprintf(buf+len, "%s%lu", sep, array[cnt]);
    len = strlen(buf);
  }
  
  char *str = init_ptr(len+1, (char*) NULL);
  strcpy(str, buf);
  delete [] buf;

  return str;
}



char *join(double *array, uint64_t size, const char *sep){
  char *buf = init_ptr(size*25, (char*) NULL);
  uint64_t len = 0;

  sprintf(buf, "%.10f", array[0]);
  len = strlen(buf);

  for(uint64_t cnt = 1; cnt < size; cnt++){
    sprintf(buf+len, "%s%.10f", sep, array[cnt]);
    len = strlen(buf);
  }

  char *str = init_ptr(len+1, (char*) NULL);
  strcpy(str, buf);
  delete [] buf;

  return str;
}


char *join(char *array, uint64_t size, const char *sep){
  char *buf = init_ptr(size*5, (char*) NULL);
  uint64_t len = 0;

  sprintf(buf, "%c", array[0]);
  len = strlen(buf);

  for(uint64_t cnt = 1; cnt < size; cnt++){
    sprintf(buf+len, "%s%c", sep, array[cnt]);
    len = strlen(buf);
  }

  char *str = init_ptr(len+1, (char*) NULL);
  strcpy(str, buf);
  delete [] buf;

  return str;
}



unsigned short int *init_ptr(uint64_t A, unsigned short int init){
  if(A < 1)
    return NULL;

  unsigned short int *ptr;
  try{
    ptr = new unsigned short int[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init;

  return ptr;
}



unsigned short int **init_ptr(uint64_t A, uint64_t B, unsigned short int init){
  if(A < 1)
    error(__FUNCTION__, "invalid size of array!");

  unsigned short int **ptr;
  try{
    ptr = new unsigned short int*[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init_ptr(B, init);

  return ptr;
}


int *init_ptr(uint64_t A, int init){
  if(A < 1)
    return NULL;

  int *ptr;
  try{
    ptr = new int[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init;

  return ptr;
}



int **init_ptr(uint64_t A, uint64_t B, int init){
  if(A < 1)
    error(__FUNCTION__, "invalid size of array!");

  int **ptr;
  try{
    ptr = new int*[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init_ptr(B, init);

  return ptr;
}



uint64_t *init_ptr(uint64_t A, uint64_t init){
  if(A < 1)
    return NULL;

  uint64_t *ptr;
  try{
    ptr = new uint64_t[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init;

  return ptr;
}



uint64_t **init_ptr(uint64_t A, uint64_t B, uint64_t init){
  if(A < 1)
    error(__FUNCTION__, "invalid size of array!");

  uint64_t **ptr;
  try{
    ptr = new uint64_t*[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init_ptr(B, init);

  return ptr;
}



float *init_ptr(uint64_t A, float init){
  if(A < 1)
    return NULL;

  float *ptr;
  try{
    ptr = new float[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init;

  return ptr;
}



double *init_ptr(uint64_t A, double init){
  if(A < 1)
    return NULL;

  double *ptr;
  try{
    ptr = new double[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init;

  return ptr;
}



double **init_ptr(uint64_t A, uint64_t B, double init){
  if(A < 1)
    error(__FUNCTION__, "invalid size of array!");

  double **ptr;
  try{
    ptr = new double*[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init_ptr(B, init);

  return ptr;
}



double ***init_ptr(uint64_t A, uint64_t B, uint64_t C, double init){
  if(A < 1)
    error(__FUNCTION__, "invalid size of array!");

  double ***ptr;
  try{
    ptr = new double**[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init_ptr(B, C, init);

  return ptr;
}



double ****init_ptr(uint64_t A, uint64_t B, uint64_t C, uint64_t D, double init){
  if(A < 1)
    error(__FUNCTION__, "invalid size of array!");

  double ****ptr;
  try{
    ptr = new double***[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init_ptr(B, C, D, init);

  return ptr;
}



// Concat str but duplicates the first one
char *strdcat(char *str1, const char *str2){
  uint64_t length = strlen(str1) + strlen(str2) + 1;

  char *str = init_ptr(length, str1);
  strcat(str, str2);

  return(str);
}



char *init_ptr(uint64_t A, const char *init){
  if(A < 1)
    return NULL;

  char *ptr;
  try{
    ptr = new char[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }
  memset(ptr, '\0', A*sizeof(char));

  if(init != NULL && strlen(init) > 0)
    strcpy(ptr, init);

  return ptr;
}



char **init_ptr(uint64_t A, uint64_t B, const char *init){
  if(A < 1)
    error(__FUNCTION__, "invalid size of array!");

  char **ptr;
  try{
    ptr = new char*[A];
  }catch (std::bad_alloc&){
    error(__FUNCTION__, "cannot allocate more memory!");
  }

  // Search for special characters
  char *pinit = init_ptr(BUFF_LEN, init);
  char *pch = strchr(pinit, '#');

  for(uint64_t a = 0; a < A; a++){
    if(pch)
      sprintf(pch, "%lu", a);
    ptr[a] = init_ptr(B, pinit);
  }
  delete [] pinit;
  delete [] pch;

  return ptr;
}



void free_ptr(void *ptr){
  delete [] (char*)ptr;
}



void free_ptr(void **ptr, uint64_t A){
  for(uint64_t a = 0; a < A; a++)
    free_ptr(ptr[a]);

  free_ptr(ptr);
}



void free_ptr(void ***ptr, uint64_t A, uint64_t B){
  for(uint64_t a = 0; a < A; a++)
    free_ptr(ptr[a], B);

  free_ptr(ptr);
}



void free_ptr(void ****ptr, uint64_t A, uint64_t B, uint64_t C){
  for(uint64_t a = 0; a < A; a++)
    free_ptr(ptr[a], B, C);

  free_ptr(ptr);
}



void cpy(void *dest, void *orig, uint64_t A, uint64_t size){
  memcpy(dest, orig, A * size);
}



void cpy(void *dest, void *orig, uint64_t A, uint64_t B, uint64_t size){
  for(uint64_t a = 0; a < A; a++)
    cpy( ((char**)dest)[a], ((char**)orig)[a], B, size);
}



void cpy(void *dest, void *orig, uint64_t A, uint64_t B, uint64_t C, uint64_t size){
  for(uint64_t a = 0; a < A; a++)
    cpy( ((char***)dest)[a], ((char***)orig)[a], B, C, size);
}



void cpy(void *dest, void *orig, uint64_t A, uint64_t B, uint64_t C, uint64_t D, uint64_t size){
  for(uint64_t a = 0; a < A; a++)
    cpy( ((char***)dest)[a], ((char***)orig)[a], B, C, D, size);
}




///////////////////////////////////
// Population Genetics Functions //
///////////////////////////////////

// Test for missing data in genotype
// => all genos are equal
bool miss_data(double *geno){
  if(abs(geno[0] - geno[1]) < EPSILON && 
     abs(geno[1] - geno[2]) < EPSILON)
    return true;
  else
    return false;
}



/* Function to call genotypes from post probs
   NOTE: probs must be normalized, that is, sum (or logsum) to 1!
      geno              : array with probs
      n_geno            : number of genotypes in array
      log_scale         : are the probs in log scale?
      N_prob_thresh     : minimum prob to use data.
                          If max prob is lower, data set as missing
      call_prob_thresh  : minimum prob to call a genotype. 
                          If max prob is lower, leave geno as probs
      miss_data         : how the missing data is handled
          0 = missing data (all genot with equal prob)
	  1 = sample random genotype
	  2 = call the highest prob geno (since missing data, probably major/major)
*/
void call_geno(double *geno, int n_geno, bool log_scale, double N_prob_thresh, double call_prob_thresh, int miss_data){
  if(N_prob_thresh > call_prob_thresh)
    error(__FUNCTION__, "missing data threshold must be smaller than calling genotype threshold!");

  int max_pos = array_max_pos(geno, n_geno);
  int min_pos = array_min_pos(geno, n_geno);
  double max_pp = (log_scale ? exp(geno[max_pos]) : geno[max_pos]);

  // If missing data
  if(geno[min_pos] == geno[max_pos]){
    if(miss_data == 0)
      max_pp = -1;
    else if(miss_data == 1)
      max_pos = rand() % 3;
  }


  if(max_pp < N_prob_thresh)
    for (int g = 0; g < n_geno; g++)
      geno[g] = (log_scale ? log((double) 1/n_geno) : (double) 1/n_geno);


  if(max_pp >= call_prob_thresh){
    for (int g = 0; g < n_geno; g++)
      geno[g] = (log_scale ? -INF : 0);

    geno[max_pos] = (log_scale ? log(1) : 1);
  }
}



// Calculate posterio probabilities (PP) from GLs and prior
// GL in log-space by default, but can be in normal-space if flag set
// prior and PP always given log-space
void post_prob(double *pp, double *lkl, double *prior, uint64_t n_geno){
  for(uint64_t cnt = 0; cnt < n_geno; cnt++){
    pp[cnt] = lkl[cnt];

    if(prior != NULL)
      pp[cnt] += prior[cnt];
  }

  double norm = logsum(pp, n_geno);

  for(uint64_t cnt = 0; cnt < n_geno; cnt++)
    pp[cnt] -= norm;
}



// Calculate HWE genotype frequencies
// MAF and F in normal-space
// Genotype frequencies in log-space
void calc_HWE(double *genot_freq, double freq, double F, bool log_scale){
  genot_freq[0] = pow(1-freq,2)   +   (1-freq)*freq*F;
  genot_freq[1] = 2*(1-freq)*freq - 2*(1-freq)*freq*F;
  genot_freq[2] = pow(freq,2)     +   (1-freq)*freq*F;

  if(log_scale)
    conv_space(genot_freq, N_GENO, log);

  /* Added to avoid impossible cases (like HET on an IBD position). This way, 
     the prior for an HET is not 0 and these cases can still be calculated. 
     We could set the PP to missing (e.g. 0.3,0.3,0.3) but the IBD status is being 
     optimized so it is probably better to keep the information and give preference to the GL.
  */
  if(F == 1){
    if(log_scale)
      genot_freq[1] = -INF;
    else
      genot_freq[1] = 1/INF;
  }
}



// Estimate site MAF from (logscaled) GL through an EM
// If indF == NULL, assumes uniform prior
// Else (0 < indF < 1), assumes HWE with specified level of inbreeding
double est_maf(uint64_t n_ind, double **pdg, double F){
  double maf;
  double *indF = init_ptr(n_ind, F);

  maf = est_maf(n_ind, pdg, indF);

  free_ptr((void*) indF);
  return maf;
}

double est_maf(uint64_t n_ind, double **pdg, double *indF){
  int iters = 0;
  double num = 0; // Expected number minor alleles
  double den = 0; // Expected total number of alleles
  double F, prev_freq, freq = 0.01;
  double prior[N_GENO], pp[N_GENO];

  do{
    prev_freq = freq;

    for(uint64_t i = 0; i < n_ind; i++){
      if(indF == NULL){
	F = 0;
	post_prob(pp, pdg[i], NULL, N_GENO);
      }else{
	F = indF[i];
	if(F < -1)
	  error(__FUNCTION__, "indF must be between 0 < F < 1!");
	calc_HWE(prior, freq, F);
	post_prob(pp, pdg[i], prior, N_GENO);
      }
      conv_space(pp, N_GENO, exp);

      num += pp[1] + pp[2]*(2-F);
      den += 2*pp[1] + (pp[0]+pp[2])*(2-F);

      //printf("Ind: %lu; num: %f; den: %f; pp: %f %f %f; IBD: %f\n", i, num, den, pp[0], pp[1], pp[2], F);
    }

    freq = num/den;
  }while( abs(prev_freq - freq) > EPSILON && iters++ < 100 );

  return freq;
}




// Adapted from BCFTOOLS:
// https://github.com/lh3/samtools/blob/6bbe1609e10f27796e5bf29ac3207bb2e35ceac8/bcftools/em.c#L266-L310
// Calculates LD directly from GLs by getting ML estimates (through an EM) of the four haplotype frequencies
void bcf_pair_LD (double LD[3], double **s1, double **s2, double maf1, double maf2, uint64_t n_ind, bool log_scale)
{
  // Estimate haplotype frequencies
  double hap_freq[4];
  haplo_freq(hap_freq, s1, s2, maf1, maf2, n_ind, log_scale);

  // Allele frequencies
  double maf[2];
  maf[0] = 1 - (hap_freq[0] + hap_freq[1]);
  maf[1] = 1 - (hap_freq[0] + hap_freq[2]);

  // calculate r^2
  double r2, D, Dp;
  D = hap_freq[0] * hap_freq[3] - hap_freq[1] * hap_freq[2]; // P_BA * P_ba - P_Ba * P_bA
  // or
  //D = hap_freq[0] - (1-maf[0]) * (1-maf[1]);               // P_BA - P_B * P_A
  Dp = D / (D < 0 ? min(maf[0]*maf[1], (1-maf[0])*(1-maf[1])) : min(maf[0]*(1-maf[1]), (1-maf[0])*maf[1]) );
  r2 = pow(D / sqrt(maf[0] * maf[1] * (1-maf[0]) * (1-maf[1])), 2);

  /*
  if(isnan(r2))
    r2 = -1.0;
  */

  LD[0] = D;
  LD[1] = Dp;
  LD[2] = r2;
}



// Get ML estimate of haplotype frequencies
int haplo_freq(double hap_freq[4], double **gl1, double **gl2, double maf1, double maf2, uint64_t n_ind, bool log_scale){
  double hap_freq_last[4];

  if(maf1 < 0 || maf1 > 1 || maf2 < 0 || maf2 > 1)
    return -1;

  // Initialize haplotype frequencies
  hap_freq[0] = (1 - maf1) * (1 - maf2); // P_BA
  hap_freq[1] = (1 - maf1) * maf2;       // P_Ba
  hap_freq[2] = maf1 * (1 - maf2);       // P_bA
  hap_freq[3] = maf1 * maf2;             // P_ba

  // iteration
  for(uint64_t i = 0; i < ITER_MAX; i++) {
    double eps = 0;
    memcpy(hap_freq_last, hap_freq, 4 * sizeof(double));
    if(log_scale)
      pair_freq_iter_log(hap_freq, gl1, gl2, n_ind);
    else
      pair_freq_iter(hap_freq, gl1, gl2, n_ind);
      
    for (uint64_t j = 0; j < 4; j++) {
      double x = fabs(hap_freq[j] - hap_freq_last[j]);
      if (x > eps) eps = x;
    }

    if(eps < EPSILON)
      break;
  }

  return 0;
}



// Function to estimate haplotype frequencies from GLs (all in normal scale)
#define _G1(h, k) ((h>>1&1) + (k>>1&1))
#define _G2(h, k) ((h&1) + (k&1))

int pair_freq_iter(double f[4], double **s1, double **s2, uint64_t n)
{
  double ff[4];
  int k, h;

  memset(ff, 0, 4 * sizeof(double));

  for (uint64_t i = 0; i < n; ++i) {
    double *p[2], sum, tmp;
    p[0] = s1[i];
    p[1] = s2[i];
    sum = 0;
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

  // Calculate frequency
  for (k = 0; k < 4; ++k)
    f[k] = ff[k] / (2 * n);

  // Normalize
  for (k = 0; k < 4; k++)
    f[k] /= f[0] + f[1] + f[2] + f[3];

  return 0;
}


int pair_freq_iter_log(double f[4], double **s1, double **s2, uint64_t n)
{
  double ff[4];
  int k, h;

  // Convert haplot frequencies to log-scale
  conv_space(f, 4, log);
  for (k = 0; k < 4; ++k)
    ff[k] = -INFINITY;

  for (uint64_t i = 0; i < n; ++i) {
    double *p[2], sum, tmp;
    p[0] = s1[i];
    p[1] = s2[i];
    sum = -INFINITY;
    for (k = 0; k < 4; ++k)
      for (h = 0; h < 4; ++h)
	sum = logsum(sum, f[k] + f[h] + p[0][_G1(k,h)] + p[1][_G2(k,h)]);
    for (k = 0; k < 4; ++k) {
      tmp = logsum(f[0] + logsum(p[0][_G1(0,k)] + p[1][_G2(0,k)], p[0][_G1(k,0)] + p[1][_G2(k,0)]),
		   f[1] + logsum(p[0][_G1(1,k)] + p[1][_G2(1,k)], p[0][_G1(k,1)] + p[1][_G2(k,1)]),
		   f[2] + logsum(p[0][_G1(2,k)] + p[1][_G2(2,k)], p[0][_G1(k,2)] + p[1][_G2(k,2)]),
		   f[3] + logsum(p[0][_G1(3,k)] + p[1][_G2(3,k)], p[0][_G1(k,3)] + p[1][_G2(k,3)]));
      ff[k] = logsum(ff[k], f[k] + tmp - sum);
    }
  }

  // Calculate frequency
  for (k = 0; k < 4; ++k)
    f[k] = exp(ff[k]) / (2 * n);

  // Normalize
  for (k = 0; k < 4; k++)
    f[k] /= f[0] + f[1] + f[2] + f[3];

  return 0;
}
