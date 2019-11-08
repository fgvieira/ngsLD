#pragma once

#include "gen_func.hpp"

// Function prototypes
double*** read_geno(char*, bool, bool, bool*, uint64_t, uint64_t);
uint64_t read_split(char *, uint64_t, uint64_t, char***, const char* = "\t");
double* read_pos(char*, uint64_t);
