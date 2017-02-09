#include "HMM.hpp"



// Forward function
double forward(double **Fw, double *q_prop, double alpha, double **e_prob, double *pos_dist, uint64_t length, int n_states){
  double *tmp = init_ptr(n_states, 0.0);
  // Initialize Fw
  for(int k = 0; k < n_states; k++)
    Fw[0][k] = log(q_prop[k]);

  // Switch from state k -> l
  for (uint64_t s = 1; s <= length; s++)
    for(int l = 0; l < n_states; l++){
      for(int k = 0; k < n_states; k++){
        tmp[k] = Fw[s-1][k] + calc_trans(k,l,q_prop[l],alpha,pos_dist[s]);

        if(isnan(tmp[k])){
          printf("site: %lu\tp_state: %d\tn_state: %d\tFw: %.15f\tq_prop: %.15f\talpha: %.15f\tdist: %f\ttrans: %f\temission: %f\n", s, k, l, Fw[s-1][k], q_prop[l], alpha, pos_dist[s], calc_trans(k,l,q_prop[l],alpha,pos_dist[s]), e_prob[s][l]);
          error(__FUNCTION__, "invalid Lkl found!");
        }
      }
      Fw[s][l] = logsum(tmp, n_states) + e_prob[s][l];
    }

  free_ptr((void*) tmp);
  return logsum(Fw[length], n_states);
}



// Backward function
double backward(double **Bw, double *q_prop, double alpha, double **e_prob, double *pos_dist, uint64_t length, int n_states){
  double *tmp = init_ptr(n_states, 0.0);
  // Initialize Bw
  for(int k = 0; k < n_states; k++)
    Bw[length][k] = log(1);

  // Switch from state k -> l
  for (uint64_t s = length; s > 0; s--){
    for(int k = 0; k < n_states; k++){
      for(int l = 0; l < n_states; l++){
	tmp[l] = calc_trans(k,l,q_prop[l],alpha,pos_dist[s]) + e_prob[s][l] + Bw[s][l];

        if(isnan(tmp[l])){
          printf("site: %lu\tp_state: %d\tn_state: %d\tBw-1: %.10f\tBw: %.10f\tq_prop: %.10f\talpha: %.10f\tdist: %f\ttrans: %f\temission: %f\n", s, k, l, tmp[l], Bw[s][l], q_prop[l], alpha, pos_dist[s], calc_trans(k,l,q_prop[l],alpha,pos_dist[s]), e_prob[s][l]);
          error(__FUNCTION__, "invalid Lkl found!");
        }
      }
      Bw[s-1][k] = logsum(tmp, n_states);
    }
  }

  // Finalize Bw
  for(int k = 0; k < n_states; k++)
    Bw[0][k] += log(q_prop[k]);

  free_ptr((void*) tmp);
  return logsum(Bw[0], n_states);
}

double viterbi_NEW(double **Vi, double *q_prop, double alpha, double **e_prob, char *path, double *pos_dist, uint64_t length, int n_states){
  int max_pos;
  int **Vi_path = init_ptr(length+1, n_states, -1);
  double *tmp = init_ptr(n_states, 0.0);

  // Initialize Vi
  for(int k = 0; k < n_states; k++)
    Vi[0][k] = log(q_prop[k]);

  for (uint64_t s = 1; s <= length; s++)
    for(int l = 0; l < n_states; l++){
      for(int k = 0; k < n_states; k++){
	tmp[k] = Vi[s-1][k] + calc_trans(k,l,q_prop[l],alpha,pos_dist[s]);

        if(isnan(tmp[k])){
          printf("site: %lu\tp_state: %d\tn_state: %d\tFw: %.15f\tq_prop: %.15f\talpha: %.15f\tdist: %f\ttrans: %f\temission: %f\n", s, k, l, Vi[s-1][k], q_prop[l], alpha, pos_dist[s], calc_trans(k,l,q_prop[l],alpha,pos_dist[s]), e_prob[s][l]);
          error(__FUNCTION__, "invalid Lkl found!");
        }
      }
      max_pos = array_max_pos(tmp, n_states);
      Vi[s][l] = tmp[max_pos] + e_prob[s][l];
      Vi_path[s][l] = max_pos;
    }


  // Trace back the Viterbi path
  path[length] = array_max_pos(Vi[length], n_states);
  for (uint64_t s = length; s > 0; s--)
    path[s-1] = Vi_path[s][(int) path[s]];

  free_ptr((void**) Vi_path, length);
  free_ptr((void*) tmp);
  return tmp[(int) path[length]];
}


double viterbi(double **Vi, double *q_prop, double alpha, double **e_prob, char *path, double *pos_dist, uint64_t length, int n_states){
  double *Vi_prob = init_ptr(n_states, 0.0);
  // Initialize Vi
  for(int k = 0; k < n_states; k++)
    Vi_prob[k] = log(q_prop[k]);

  for (uint64_t s = 1; s <= length; s++){
    for(int l = 0; l < n_states; l++){
      double vmax = -INF;
      uint64_t k_vmax = 0;

      for(int k = 0; k < n_states; k++){
        double pval = Vi_prob[k] + calc_trans(k,l,q_prop[l],alpha,pos_dist[s]);
	if( vmax < pval ) { vmax = pval; k_vmax = k; }
      }

      Vi[s][l] = k_vmax;
      Vi_prob[l] = vmax + e_prob[s][l];
    }
  }

  // Trace back the Viterbi path
  path[length] = array_max_pos(Vi_prob, n_states);
  for (uint64_t s = length; s > 0; s--)
    path[s-1] = Vi[s][(int) path[s]];

  return Vi_prob[(int) path[length]];
}



// Calculates transition probabilities between states k and l, depending on distance, inbreeding and transition rate
double calc_trans(char k, char l, double q_prop, double alpha, double pos_dist){
  double trans = 0;
  double coanc_change = exp(-alpha*pos_dist);
  
  trans = (1-coanc_change) * q_prop;
  if(k == l)
    trans += coanc_change;

  return log(trans);
}



// Calculates emission probabilities depending on the MAF and inbreeding level (state).
double calc_emission(double gl[3], double maf, uint64_t F){
  if(maf < 0 || maf > 1)
    error(__FUNCTION__, "invalid MAF!");
  if(F < 0 || F > 1)
    error(__FUNCTION__, "invalid F!");

  double geno[3];
  calc_HWE(geno, maf, F, true);

  return logsum(gl[0]+geno[0], gl[1]+geno[1], gl[2]+geno[2]);
}

double calc_emission(double gl[3], double maf1, double maf2){
  if(maf1 < 0 || maf1 > 1)
    error(__FUNCTION__, "invalid MAF1!");
  if(maf2 < 0 || maf2 > 1)
    error(__FUNCTION__, "invalid MAF2!");

  double geno[3];
  geno[0] = log((1 - maf1) * (1 - maf2));
  geno[1] = log((1-maf1) * maf2 + maf1 * (1-maf2));
  geno[2] = log(maf1 * maf2);

  return logsum(gl[0]+geno[0], gl[1]+geno[1], gl[2]+geno[2]);
}



// Calculates emission probabilities depending on the MAF and inbreeding level (state).
// site _p is previous site
// site _c is current site
double calc_emissionLD(double hap_freq[4], double *gl_p, double *gl_c, double maf_p, double maf_c,  uint64_t F){
  if(maf_p < 0 || maf_p > 1 || maf_c < 0 || maf_c > 1)
    error(__FUNCTION__, "invalid MAF!");

  double sum, P_Gc, P_Gp[N_GENO];
  double s_p[N_GENO], s_c[N_GENO];
  for(uint64_t g = 0; g < N_GENO; g++){
    s_p[g] = exp(gl_p[g]);
    s_c[g] = exp(gl_c[g]);
  }

  sum = 0;
  if(0){
    for(uint64_t g_c = 0; g_c < N_GENO; g_c++){
      calc_HWE(P_Gp, maf_p, F);
      for(uint64_t g_p = 0; g_p < N_GENO; g_p++)
	P_Gc += joint_geno_prob(hap_freq, g_p, F, g_c, F) * P_Gp[g_p];

      sum += s_c[g_c] * P_Gc;
    }

    return log(sum);
  }else{
    // results are calculated with current algorithm
    // results _2 are calculated just taking the haplotype frequency into account (removing s_p[g_p] * s_c[g_c])

    for(uint64_t g_c = 0; g_c < N_GENO; g_c++)
      for(uint64_t g_p = 0; g_p < N_GENO; g_p++)
	sum += joint_geno_prob(hap_freq, g_p, F, g_c, F) * s_p[g_p] * s_c[g_c];

    return log(sum) - calc_emission(gl_p, maf_p, F);
  }
}



double joint_geno_prob(double hap_freq[4], uint64_t g_p, uint64_t F_p, uint64_t g_c, uint64_t F_c){
  if(F_p != F_c)
    error(__FUNCTION__, "so far only cases where prev and curr positions have same inbreeding level!");

  if(g_p == 0 && g_c == 0){
    return (F_c == 0 ? pow(hap_freq[0], 2) : hap_freq[0]);
  }else if(g_p == 0 && g_c == 1){
    return (F_c == 0 ? 2*hap_freq[0]*hap_freq[1] : 0);
  }else if(g_p == 0 && g_c == 2){
    return (F_c == 0 ? pow(hap_freq[1], 2): hap_freq[1]);
  }else if(g_p == 1&& g_c == 0){
    return (F_c == 0 ? 2*hap_freq[0]*hap_freq[2] : 0);
  }else if(g_p == 1&& g_c == 1){
    return (F_c == 0 ? 2*(hap_freq[0]*hap_freq[3]+hap_freq[1]*hap_freq[2]) : 0);
  }else if(g_p == 1 && g_c == 2){
    return (F_c == 0 ? 2*hap_freq[1]*hap_freq[3] : 0);
  }else if(g_p == 2 && g_c == 0){
    return (F_c == 0 ? pow(hap_freq[2], 2) : hap_freq[2]);
  }else if(g_p == 2 && g_c == 1){
    return (F_c == 0 ? 2*hap_freq[2]*hap_freq[3] : 0);
  }else if(g_p == 2 && g_c == 2){
    return (F_c == 0 ? pow(hap_freq[3], 2) : hap_freq[3]);
  }

  return -1;
}
