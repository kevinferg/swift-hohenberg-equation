#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*  _________________________________  */
/* |__________  Constants  __________| */

#define RES                     32     // Number of rows (and columns)
static const float scale      = 32.0f; // Number of cycles across domain
static const float dt         = 0.2f;  // Time step
static const int   num_steps  = 250;   // Number of iterations
static const float epsilon    = 1.0f;  // Linear growth coefficient
static const float wavenum    = 1.0f;  // Wave number

/*  _________________________________  */


static float u[RES*RES] = {0};

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

void print_array(FILE* stream, float* array, int rows, int cols, float vmin, float vmax, char* chars) {
    int r, c, index;
    int n = strlen(chars);
    if (n <= 0) return;
    float range = vmax - vmin;
    float h = range / (float) n;
    float val;
    char p;
    for (r = 0; r < rows; r++) {
        for (c = 0; c < cols; c++) {
            val = array[cols*r + c];
            index = (int) ((val - vmin) / range * n);
            index = index < 0? 0: index > n-1? n-1: index;
            p = chars[index];
            // printf("%f %d\n", (val - minval) / range * n, (int) ((val - minval) / range * n));
            fprintf(stream,"%c", p);
        }
        fprintf(stream, "\n");
    }
    return;
}

int main(int argc, char** argv) {
    random_normal_array(u, RES*RES, 0, 1);
    print_array(stdout, u, RES, RES, -3, 3, " .oO@");
    return 0;
}