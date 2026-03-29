#ifndef SH_EQN_H
#define SH_EQN_H

#include "stdint.h"

typedef struct SHOptions {
    float scale;         // Waviness
    float dt;            // Time step
    float epsilon;       // Linear growth coefficient
    float wavenum;       // Wave number
    float init_stdev;    // Stdev of initial noise field
    int32_t num_steps;   // Number of iterations
} SHOptions;

extern SHOptions default_sh_options;

// WARNING: This mallocs the result, so be sure to free it!
// Also, if you supply res other than a power of 2, it gets rounded to one of them.
char* generate_sh_field(SHOptions* options, int32_t res, uint32_t seed, int32_t charwidth, char* charmap);

#endif