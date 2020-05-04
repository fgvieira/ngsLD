#pragma once

#include "gen_func.hpp"

// Function prototypes
double*** read_geno(char*, bool, bool, bool*, uint64_t, uint64_t);
int read_split(char *, char***, uint64_t*, uint64_t*, const char* = "\t");
double* read_dist(char*, uint64_t);
