#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>
typedef float complex complex32_t;

#include "sh-eqn.h"

/**********************************************************************/
/*                             Constants                              */
/**********************************************************************/
#define MAX_RES 512
SHOptions default_sh_options = {
    .scale = 1.0f,
    .dt = 0.2f,
    .num_steps = 200,
    .epsilon = 1.0f,
    .wavenum = 1.0f,
    .init_stdev = 0.1f,
};
static uint32_t rng_seed = 321;
/*                                                                    */
/**********************************************************************/

// 2D float array utility functions
static void random_normal_array(float* arr, int N, float mean, float stdev);
void print_array(FILE* stream, float* array, int rows, int cols, float vmin, float vmax, char* chars, int width);

// FFT functions
void fft_1d_c32(complex32_t *x, int n);
void fft2_real(float *input_real, complex32_t *output_freq, complex32_t *buf, int n);
void ifft2_complex(complex32_t *input_freq, float *output_real, complex32_t *buf, int n);

// Swift Hohenberg solver
int solve_swift_hohenberg(float* u, int res, SHOptions* options);

/********************************************************************************/
/*                     2D float array utility functions                         */
/********************************************************************************/

static inline uint16_t rand15(void) {
    rng_seed = rng_seed * 747796405u + 2891336453u;
    return rng_seed >> 17;  // return 15-bit result
}

static inline float frandn(void) {
    static const int16_t samples[5][8] = {
        {-18666, 1465,  9124, -2456,  7113,   2465,  2440,  -4190},
        { -1690, 1656, -6392,  3616, -7014,   8614,  5799,  -1058},
        { -3021, 4302, -3768, -7831, 11343,  13204, -1676,  -5054},
        {  1012, 6001,  5497, 14892,   487,    967,   -39, -13798},
        { -7858, 3286,  1206, -5579, 14569, -11125, -9127,  -8716}
    }; // ^ Fine-tuned samples from ~N(mu=0, sig=2^14/sqrt(5))
    uint16_t r = rand15();          // 15 random bits
    int32_t x;                      // Add 1 entry from each row
    x  = samples[0][r&7];  r >>= 3;
    x += samples[1][r&7];  r >>= 3;
    x += samples[2][r&7];  r >>= 3;
    x += samples[3][r&7];  r >>= 3;
    x += samples[4][r&7];           // Now, x~N(mu=0, sig=2^14)
    return (float) x * (1.0f / 16384.0f);
}

void random_normal_array(float* arr, int N, float mean, float stdev) {
    uint32_t i;
    for (i = 0; i < N; i++) arr[i] = frandn() * stdev + mean;
}

void print_array(FILE* stream, float* array, int rows, int cols, float vmin, float vmax, char* chars, int width) {
    int r, c, index, j, n;
    if ( (n = strlen(chars)) <= 0) return;
    if (vmin == vmax) {
        vmin = array[0];
        vmax = array[1];
        for (j = 1; j < n; j++) {
            if (array[j] < vmin) vmin = array[j];
            if (array[j] > vmax) vmax = array[j];
        }
    }
    float range = vmax - vmin;
    for (r = 0; r < rows; r++) {
        for (c = 0; c < cols; c++) {
            index = (int) ((array[cols*r + c] - vmin) / range * n);
            index = index < 0? 0: index > n-1? n-1: index;
            for (j = 0; j < width; j++) 
                fprintf(stream,"%c", chars[index]);
        }
        fprintf(stream, "\n");
    }
}

char* sprint_array(float* array, int rows, int cols, float vmin, float vmax, char* chars, int width) {
    int r, c, index, j, n;
    if ((n = strlen(chars)) <= 0) return NULL;

    if (vmin == vmax) {
        vmin = array[0];
        vmax = array[1];
        for (j = 1; j < n; j++) {
            if (array[j] < vmin) vmin = array[j];
            if (array[j] > vmax) vmax = array[j];
        }
    }

    float range = vmax - vmin;

    int line_len = cols * width + 1;
    int total_len = rows * line_len + 1;
    char* out = malloc(total_len);
    char* p = out;

    for (r = 0; r < rows; r++) {
        for (c = 0; c < cols; c++) {
            index = (int)((array[cols*r + c] - vmin) / range * n);
            index = index < 0 ? 0 : index > n-1 ? n-1 : index;
            for (j = 0; j < width; j++)
                *p++ = chars[index];
        }
        *p++ = '\n';
    }

    *p = '\0';
    return out;
}

/********************************************************************************/
/*                                 FFT Functions                                */
/********************************************************************************/

void fft_1d_c32(complex32_t *x, int n) {
    // Note: n must be a power of 2.

    int index, bit_reversed_index;
    for (index = 1, bit_reversed_index = 0; index < n; index++) {
        int bit_mask = n >> 1;
        while (bit_reversed_index & bit_mask) {
            bit_reversed_index &= ~bit_mask;
            bit_mask >>= 1;
        }

        bit_reversed_index |= bit_mask;
        if (index < bit_reversed_index) {
            complex32_t temp = x[index];
            x[index] = x[bit_reversed_index];
            x[bit_reversed_index] = temp;
        }
    }

    int stage_size, half_size, k;
    complex32_t wlen, w, rotated_value;
    for (stage_size = 2; stage_size <= n; stage_size <<= 1) {
    // Iterate over FFT stages, sub-FFT sizes doubling each time
        half_size = stage_size >> 1;
        wlen = cexpf(-2.0f * I * M_PI / stage_size);

        for (int block_start = 0; block_start < n; block_start += stage_size) {
        // Iterate over independent butterfly blocks within this stage
            w = 1.0f;
            for (k = 0; k < half_size; k++) {
            // Individual butterfly operations within a block
                int lower_index = block_start + k;
                int upper_index = lower_index + half_size;
                rotated_value = w * x[upper_index];
                x[upper_index] = x[lower_index] - rotated_value;
                x[lower_index] = x[lower_index] + rotated_value;
                w *= wlen;
            }
        }
    }
}

void fft2_real(float *input_real, complex32_t *output_freq, complex32_t *buf, int n) {
    int x, y;
    for (y = 0; y < n; y++) { // Row FFTs
        int row_offset = y * n;
        for (x = 0; x < n; x++) buf[x] = input_real[row_offset + x];
        fft_1d_c32(buf, n);
        for (x = 0; x < n; x++) output_freq[row_offset + x] = buf[x];
    }
    for (x = 0; x < n; x++) { // Column FFTs
        for (y = 0; y < n; y++) buf[y] = output_freq[y * n + x];
        fft_1d_c32(buf, n);
        for (y = 0; y < n; y++) output_freq[y * n + x] = buf[y];
    }
}

void ifft2_complex(complex32_t *input_freq, float *output_real, complex32_t *buf, int n) {
    int x, y;
    for (y = 0; y < n; y++) { // Row inverse FFTs
        int row_offset = y * n;
        for (x = 0; x < n; x++) buf[x] = conjf(input_freq[row_offset + x]);
        fft_1d_c32(buf, n);
        for (x = 0; x < n; x++) input_freq[row_offset + x] = conjf(buf[x]);
    }

    for (x = 0; x < n; x++) { // Column inverse FFTs
        for (y = 0; y < n; y++) buf[y] = conjf(input_freq[y * n + x]);
        fft_1d_c32(buf, n);
        for (y = 0; y < n; y++) output_real[y * n + x] = crealf(buf[y]) / (n * n);
    }
}

/********************************************************************************/
/*                          Swift-Hohenberg Solver                              */
/********************************************************************************/

int32_t get_valid_res(int32_t x) {
    if (x <= 4)
        return 4;

    x |= x >> 1;   x |= x >> 2;
    x |= x >> 4;   x |= x >> 8;
    x |= x >> 16;  x = x - (x >> 1);

    if (x > MAX_RES) return MAX_RES;
    return x;
}

int solve_swift_hohenberg(float* u, int res, SHOptions* options) {
    /*  
        u'(t) = epsilon*u - (Del^2 + wavenum^2)^2 * u - u^3

                ____________Linear Operator__________     ___Nonlinear___
        u'(t) = [epsilon - (Del^2 + wavenum^2)^2] * u  +       -u^3
              =                 L(u)                   +       N(u)

        Euler: u(t+1) = u(t) + dt * u'()
                      = u(t) + dt*L(u(t+1))  +  dt*N(u(t))
               u(t+1) - dt*L(u(t+1)) = u(t)  +  dt*N(u(t))

        Fourier Transform.... Del^2 becomes multiplication by -k^2
        (1 - dt*(epsilon - (wavenum^2 - k^2)^2)) * FFT[u(t+1)] = FFT[u(t) + dt*N(u(t))]
        ... Precompute lin. op. array:  denom = 1 - dt*(epsilon - (wavenum^2 - k^2)^2)
        FFT[u(t+1)] = FFT[u(t) + dt*N(u(t))] / denom

        u(t+1) = iFFT{FFT[u(t) + dt*N(u(t))] / denom}
    */
    int i, x, y;
    size_t total_size = 2*res*res*sizeof(complex32_t) + (2*res*res + res)*sizeof(float);
    void *mem = malloc(total_size);
    if (!mem) return -1;

    complex32_t *uhat = mem;
    complex32_t *uhat_buf = uhat + res*res;
    float *denom = (float *)(uhat_buf + res*res);
    float *real_vals = denom + res*res;
    float *K2 = real_vals + res*res;

    float lin_op, dk = 2.0f * M_PI / (options->scale * res);
    for (i = 0; i <= res/2; i++) {
        K2[i] = i * dk;
        K2[i] *= K2[i]; // Laplacian becomes elementwise square of meshgrid in freq space
        K2[res - i] = K2[i]; // Values are mirrored + negated across middle index
        // (Skipping negation step bc it's all squared)
    }
    for (y = 0; y < res; y++) {
        for (x = 0; x < res; x++) {
            lin_op = (options->wavenum*options->wavenum - (K2[x] + K2[y]));
            lin_op = options->epsilon - lin_op*lin_op;
            denom[y*res + x] = 1.0f - options->dt*lin_op;
        }
    }

    for (i = 0; i < options->num_steps; i++) {
        for (y = 0; y < res*res; y++)
            u[y] = u[y] + options->dt * ( -u[y]*u[y]*u[y] ); // Real space operations

        fft2_real(u, uhat, uhat_buf, res);

        for (y = 0; y < res*res; y++)
            uhat[y] /= denom[y]; // Elementwise operations in freq space

        ifft2_complex(uhat, u, uhat_buf, res);
    }
    free(mem);
    return 0;
}

/********************************************************************************/
/*                         Main user-facing functions                           */
/********************************************************************************/


char* generate_sh_string(SHOptions* options, int32_t res, uint32_t seed, int32_t charwidth, char* charmap) {
    res = get_valid_res(res);
    float* u = generate_sh_field(options, res, seed);
    if (!u) return NULL;
    char* str = sprint_array(u, res, res, 0, 0, charmap, charwidth);
    free(u);
    return str;
}

float* generate_sh_field(SHOptions* options, int32_t res, uint32_t seed) {
    res = get_valid_res(res);
    if (options == NULL) {
        options = &default_sh_options;
    }
    rng_seed = seed;
    float* u;
    u = calloc(res*res, sizeof(float));
    if (!u) return NULL;
    random_normal_array(u, res*res, 0, options->init_stdev);
    solve_swift_hohenberg(u, res, options);
    return u;
}