#include <stdio.h>
#include <stdlib.h>

#include "sh-eqn.h"

int main(int argc, char** argv) {
    SHOptions options = default_sh_options;
    options.scale = 1.25f;
    // options.dt = 0.2f;
    // options.num_steps = 200;
    // options.epsilon = 1.0f;
    // options.wavenum = 1.0f;
    // options.init_stdev = 0.1f;

    char* example_solution_string;
    example_solution_string = generate_sh_field(
        &options, // Misc. parameters
        32,       // Resolution (power of 2)
        234,      // RNG seed
        2,        // Print 2 chars per cell
        " .:o0@"  // ASCII chars to print
    );
    puts(example_solution_string);

    free(example_solution_string);
    return 0;
}