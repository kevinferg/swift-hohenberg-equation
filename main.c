#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
typedef float complex complex32_t;
#ifndef PI
#define PI 3.14159265358979
#endif


/**********************************************************************/
/*                             Constants                              */
/**********************************************************************/
#define RES                     64     // Number of rows (and columns)
static const float scale      = 64.0f; // Number of cycles across domain
static const float dt         = 0.2f;  // Time step
static const int num_steps    = 200;   // Number of iterations
static const float epsilon    = 1.0f;  // Linear growth coefficient
static const float wavenum    = 1.0f;  // Wave number
static const float init_stdev = 0.1f;  // Stdev of initial noise field
static const int print_width  = 2;     // Number of chars per float
/*                                                                    */
/**********************************************************************/

// 2D float array utility functions
void random_normal_array(float* arr, int N, float mean, float stdev);
void print_array(FILE* stream, float* array, int rows, int cols, float vmin, float vmax, char* chars, int width);

// FFT functions
void fft_1d_c32(complex32_t *x, int n);
void fft2_real(float *input_real, complex32_t *output_freq, complex32_t *buf, int n);
void ifft2_complex(complex32_t *input_freq, float *output_real, complex32_t *buf, int n);

// Swift Hohenberg solver
int solve_swift_hohenberg(float* u, int res);



int main(int argc, char** argv) {
    static float u[RES*RES] = {0};
    random_normal_array(u, RES*RES, 0, init_stdev);
    solve_swift_hohenberg(u, RES);
    print_array(stdout, u, RES, RES, 0, 0, " .:o0@", print_width);
    return 0;
}

/********************************************************************************/
/*                     2D float array utility functions                         */
/********************************************************************************/

void random_normal_array(float* arr, int N, float mean, float stdev) {
    static const int16_t prob0[16] = {-9555, -5842, -3891, -2666, -1954, -1112, -713, 370, -265, 339, 1480, 2780, 2508, 3849, 6098, 8303};
    static const int16_t prob1[16] = {-8339, -5311, -4676, -2430, -1969, -1431, -101, -284, -225, 968, 1615, 2042, 2668, 4516, 6143, 9003};
    static const int16_t prob2[16] = {-8303, -5924, -3739, -2536, -2945, -1728, -1137, -1282, 339, 1408, 1334, 1334, 3020, 3988, 4477, 8564};
    static const int16_t prob3[8]  = {-6486, -3190, -1834, -801, 1022, 2116, 2882, 6897};
    int i, x;
    unsigned int r;
    for (i = 0; i < N; i++) {
        r = rand();
        x  =  prob0[r&15]; r >>= 4;
        x +=  prob1[r&15]; r >>= 4;
        x +=  prob2[r&15]; r >>= 4;
        x +=  prob3[r&7];

        arr[i] = mean + stdev * ((float) x / (float) (8192.f));
    }
    return;
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
    return;
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
        wlen = cexpf(-2.0f * I * PI / stage_size);

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
    return;
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
    return;
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

int solve_swift_hohenberg(float* u, int res) {
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

    float lin_op, dk = 2.0f * PI / scale;
    for (i = 0; i <= RES/2; i++) {
        K2[i] = i * dk;
        K2[i] *= K2[i]; // Laplacian becomes elementwise square of meshgrid in freq space
        K2[RES - i] = K2[i]; // Values are mirrored + negated across middle index
        // (Skipping negation step bc it's all squared)
    }
    for (y = 0; y < RES; y++) {
        for (x = 0; x < RES; x++) {
            lin_op = (wavenum*wavenum - (K2[x] + K2[y]));
            lin_op = epsilon - lin_op*lin_op;
            denom[y*RES + x] = 1.0f - dt*lin_op;
        }
    }

    for (i = 0; i < num_steps; i++) {
        for (y = 0; y < RES*RES; y++)
            u[y] = u[y] + dt * ( -u[y]*u[y]*u[y] ); // Real space operations

        fft2_real(u, uhat, uhat_buf, RES);

        for (y = 0; y < RES*RES; y++)
            uhat[y] /= denom[y]; // Elementwise operations in freq space

        ifft2_complex(uhat, u, uhat_buf, RES);
    }
    free(mem);
    return 0;
}