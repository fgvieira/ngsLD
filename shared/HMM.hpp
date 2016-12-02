#pragma once

#include "gen_func.hpp"



// Function prototypes
double forward(double**, double*, double, double**, double*, uint64_t, int);
double backward(double**, double*, double, double**, double*, uint64_t, int);
double viterbi(double**, double*, double, double**, char*, double*, uint64_t, int);
double calc_trans(char, char, double, double, double);
double calc_emission(double [N_GENO], double, uint64_t);
double calc_emissionLD(double [4], double*, double*, double, double, uint64_t);
double joint_geno_prob(double [4], uint64_t, uint64_t, uint64_t, uint64_t);
