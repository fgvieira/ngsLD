#include "gen_func.hpp"


bool SIG_COND;
int really_kill = 3;


void warn(const char *func, const char *msg) {
  fflush(stdout);
  fprintf(stderr, "\n[%s] WARN: %s\n", func, msg);
  fflush(stderr);
}


void error(const char *func, const char *msg) {
  fflush(stdout);
  fprintf(stderr, "\n[%s] ERROR: %s\n", func, msg);
  perror("\t");
  fflush(stderr);
  exit(-1);
}


void handler(int s) {
  if(SIG_COND)
    fprintf(stderr,"\n\"%s\" signal caught! Will try to exit nicely (no more threads are created but we will wait for the current threads to finish).\n", strsignal(s));

  if(--really_kill != 3)
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
  double errTol = 1e-5;

  if (value != value) {
    error(__FUNCTION__, "value is NaN!\n");
  } else if(value < errTol) {
    value = 0;
    if(verbose && value < 0)
      fprintf(stderr, "\nWARN: value %f < 0!\n", value);
  } else if(value > 1 - errTol) {
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
void transp_matrix(void *src, void *dest, uint64_t A, uint64_t B){
  if(A < 1 || B < 1)
    error(__FUNCTION__, "invalid size of array!");

  double** s = (double**) src;
  double** d = (double**) dest;
  
  for(uint64_t a = 0; a < A; a++)
    for(uint64_t b = 0; b < B; b++)
      d[b][a] = s[a][b];
}



double draw_rnd(gsl_rng *r, uint64_t min, uint64_t max) {
  return( min + gsl_rng_uniform(r) * (max - min) );
}



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



void conv_space(double *geno, int n_geno, double (*func)(double)) {
  for (int g = 0; g < n_geno; g++) {
    geno[g] = func(geno[g]);

    if(geno[g] == -INFINITY)
      geno[g] = -INF;
  }
}



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
double logsum2(double a, double b){
  double buf[2];
  buf[0] = a;
  buf[1] = b;
  return logsum(buf, 2);
}


// Special logsum case for size == 3
double logsum3(double a, double b, double c){
  double buf[3];
  buf[0] = a;
  buf[1] = b;
  buf[2] = c;
  return logsum(buf, 3);
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
uint64_t read_file(const char *in_file, char ***ptr, uint64_t buff_size){
  uint64_t cnt = 0;
  char buf[buff_size];
  char **tmp = NULL;
  uint64_t max_len = 0;

  // Open file
  gzFile in_file_fh = gzopen(in_file, "r");
  if(in_file_fh == NULL)
    error(__FUNCTION__, "cannot open file");

  while(!gzeof(in_file_fh)){
    buf[0] = '\0';
    // Read line from file
    gzgets(in_file_fh, buf, buff_size);
    // Remove trailing newline
    chomp(buf);
    // Check if empty
    if(strlen(buf) == 0)
      continue;
    // Alloc memory
    tmp = (char**) realloc(tmp, (cnt+1)*sizeof(char*));
    tmp[cnt] = (char*) calloc(strlen(buf)+1, sizeof(char));
    strcpy(tmp[cnt], buf);
    cnt++;
    if(strlen(buf) > max_len)
      max_len = strlen(buf);
  }

  // Copy to final array
  *ptr = init_ptr(cnt, max_len+1, (const char*) '\0');
  for(uint64_t i = 0; i < cnt; i++){
    strcpy(ptr[0][i], tmp[i]);
    free(tmp[i]);
  }
  free(tmp);
  
  gzclose(in_file_fh);
  return cnt;
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
uint64_t split(char *str, const char *sep, int **out){
  uint64_t i = strlen(str);
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


uint64_t split(char *str, const char *sep, float **out){
  uint64_t i = strlen(str);
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


uint64_t split(char *str, const char *sep, double **out){
  uint64_t i = strlen(str);
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


uint64_t split(char *str, const char *sep, char ***out){
  uint64_t i = strlen(str);
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

  sprintf(buf, "%d", array[0]);
  len = strlen(buf);

  for(uint64_t cnt = 1; cnt < size; cnt++){
    sprintf(buf+len, "%s%d", sep, array[cnt]);
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
  for(uint64_t a = 0; a < A; a++)
    ptr[a] = init_ptr(B, init);

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
