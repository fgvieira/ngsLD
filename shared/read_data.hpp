#pragma once

#include "gen_func.hpp"


const uint64_t N_GENO = 3;


double*** read_geno(char*, bool, bool, bool, uint64_t, uint64_t, bool = false);
double* read_pos(char*, uint64_t);
