# swift-hohenberg-equation
Solving the Swift-Hohenberg equation to generate stripe patterns

$$
u'(t) = \varepsilon u - (\lambda^2 + \nabla^2)^2 \cdot u - u^3
$$

## Usage

Compile (on Windows with gcc installed) with the following:
```
gcc -O2 -lm *.c -o main.exe
```


## Example

Compiling and running with these parameters:
```c
#define RES                     32     // Number of rows (and columns)
static const float scale      = 32.0f; // Number of cycles across domain
static const float dt         = 0.2f;  // Time step
static const int num_steps    = 200;   // Number of iterations
static const float epsilon    = 1.0f;  // Linear growth coefficient
static const float wavenum    = 1.0f;  // Wave number
static const float init_stdev = 0.1f;  // Stdev of initial noise field
static const int print_width  = 2;     // Number of chars per float
static uint16_t rng_seed      = 321;   // RNG seed
```

... results in the following being printed to the console:
```
  ..00@@@@oo    ..00@@@@oo..    oo@@@@00::......    ..oo@@@@oo  
    oo@@@@oo    ::@@@@00::    ::00@@00::    ::oo::..  ::00@@00..
    oo@@@@oo    ::@@@@00..    oo@@@@oo    ..oo@@00..  ..00@@@@::
    ::@@@@oo    ::@@@@oo    ..00@@@@::    ..00@@@@..    oo@@@@::
    oo@@@@oo    ::@@@@oo..  ..oo@@@@oo    ::@@@@00..  ..00@@@@::
    oo@@@@oo    ::@@@@00..    ::00@@00::::00@@@@oo    ..00@@00..
  ..00@@@@oo    ::00@@@@oo    ..::00000000@@@@oo..    oo@@@@oo  
  ::@@@@00::    ..oo@@@@00::      ::oo@@@@@@00::    ::00@@00::
..oo@@@@oo..      ::00@@@@00::      ::00@@@@oo      oo@@@@oo..
::00@@00::    ......::oo@@@@00::    ..oo@@@@::    ..00@@@@::
oo@@@@oo    ::oo::..  ..oo@@@@00..    oo@@@@::    ::@@@@00..
00@@00..  ..oo@@00::    ..00@@@@::    oo@@@@oo    ::@@@@oo..  ..
00@@00..  ..00@@@@::    ..00@@@@::    oo@@@@oo    ::@@@@oo    ..
00@@00..  ..oo@@@@oo    ..00@@@@::    oo@@@@::    ::@@@@oo    ..
00@@00::    ::@@@@00::::oo@@@@00..    oo@@@@::    ::@@@@oo    ..
oo@@@@oo    ..oo@@@@000000@@00::    ..oo@@@@::    ::@@@@00..  ..
::@@@@00..    ::00@@@@@@00oo::      ::00@@00::    ::00@@00..
..00@@@@::    ..00@@@@00::..    ..::oo@@@@@@::    ..00@@@@::
..00@@@@oo    ..00@@@@::      ::oo0000@@@@@@oo..    oo@@@@oo    
  oo@@@@oo    ..00@@00..    ::00@@@@00oo00@@00::    ::@@@@00..
..oo@@@@::    ::@@@@00..  ..oo@@@@00::..oo@@@@::    ::@@@@00..
..00@@@@::    ::@@@@oo    ..00@@@@::    ::@@@@oo    ::@@@@00..
..00@@00..    oo@@@@oo    ::@@@@00..    ::@@@@oo    ::@@@@oo
::@@@@oo    ..00@@@@::    ::@@@@oo      oo@@@@::    oo@@@@oo
oo@@@@oo    ::@@@@00..    oo@@@@oo    ..00@@00..    oo@@@@::
oo@@@@::    oo@@@@00..  ..oo@@@@::    ::@@@@00..  ..00@@00..    
00@@00..    oo@@@@oo    ..00@@00::    ::@@@@oo    ..00@@00..  ..
00@@00::    oo@@00::    ::@@@@00..    oo@@@@oo    ::@@@@00..  ..
oo@@@@::    ::oo::..  ..00@@@@::    ..00@@@@::    ::@@@@00..
::@@@@00..    ....  ..oo@@@@oo..    ::@@@@00::    ::@@@@00..
..00@@@@oo..      ..oo@@@@00::    ::00@@@@oo..    ::00@@@@oo
  ::@@@@00::    ..oo@@@@00::    ..oo@@@@00::      ..oo@@@@00..
```


## Explanation

The Swift-Hohenberg Equation is given as  
$$
u'(t) = \varepsilon u - (\lambda^2 + \nabla^2)^2 \cdot u - u^3.
$$

The right-hand side can be rewritten as a sum of a linear operator and a nonlinear operator on $u$:  
$$
u'(t) =  \quad   [\varepsilon - (\lambda^2 + \nabla^2)^2] \cdot u   \quad  +   \quad  -u^3  \quad = L(u) + N(u).
$$

Starting with an initial condition $u(0)$ given, we can integrate with Euler's method:  
$$
u(t+1) = u(t) + dt * u'(t).
$$

We make the method semi-implicit by representing $u'(t) = L(u(t+1)) + N(u(t))$, such that:  
$$
u(t+1) = u(t) + dt \cdot L(u(t+1))  +  dt \cdot N(u(t)).
$$

Rearranging, we get the u(t+1) terms on one side and the u(t) terms on the other:  
$$
u(t+1) - dt \cdot L(u(t+1)) = u(t) +  dt \cdot N(u(t)).
$$

Now, we take the Fourier Transform of both sides. The $\nabla^2\cdot$ operator becomes multiplication by -k^2 in the Fourier domain (for the 2D problem, this is a grid of frequency values). The resulting equation after applying an FFT is thus:  
$$
(1 - dt\cdot (\varepsilon - (\lambda^2 - k^2)^2)) * FFT[u(t+1)] = FFT[u(t) + dt\cdot N(u(t))].
$$

For speed, we can precompute the linear operator array $Q = 1 - dt \cdot (\varepsilon - (\lambda^2 - k^2)^2)$ across the meshgrid of frequencies $k$. This gives the following in Fourier space:  
$$
FFT[u(t+1)] = FFT[u(t) + dt \cdot N(u(t))] / Q.
$$

Thus, a full semi-implicit Euler time step update can be written as:  
$$
u(t+1) = iFFT\{FFT[u(t) + dt\cdot N(u(t))] / Q\}
$$

In practice, this means:
1. Set initial grid `u(0)` as Gaussian noise
2. Precompute meshgrid of `k` frequencies and the denominator term `Q = 1 - dt*(epsilon - (wavenum^2 - k^2)^2)`
3. Starting with current grid `u`, evaluate `u + dt * N(u)`
4. Compute the 2D FFT
5. Divide elementwise by the precomputed denominator term `Q`
6. Compute the 2D inverse FFT
7. Repeat steps 3-6 for desired number of time steps