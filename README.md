# swift-hohenberg-equation
Solving the Swift-Hohenberg equation to generate stripe patterns

$$
u'(t) = \varepsilon u - (\lambda^2 + \nabla^2)^2 \cdot u - u^3
$$

## Usage

Compile using `make` or the following command:
```
gcc -O2 -lm *.c -o main.exe
```

Then running `./main.exe` results in the following being printed to the console:  


```
@@00..  ::@@@@..    oo@@@@00..  ..oo@@@@@@@@::    00@@@@::  ..00
@@00    ::@@@@::      oo@@@@00..    ::@@@@00..  ..@@@@00    ::@@
@@00..  ..@@@@@@::    ..oo@@@@00..    oo@@@@::  ..@@@@00    ::@@
@@@@oo    ::@@@@@@00..    oo@@@@00..  ..@@@@oo  ..00@@@@::    oo
00@@@@::    ::@@@@@@00::    oo@@@@oo    00@@00..  ::@@@@00..  ..
..@@@@@@..    ::00@@@@00..  ..00@@00..  oo@@@@::    oo@@@@00..  
  ::@@@@00::      ::00@@00    ::@@@@::  ..@@@@00    ..00@@@@oo  
    oo@@@@@@oo      ::@@@@::  ..@@@@oo    00@@@@oo    ..00@@@@::
..  ..oo@@@@@@00....oo@@@@oo    00@@00..  ::@@@@@@oo    ..00@@@@
oo      ::@@@@@@@@@@@@@@@@::  ..@@@@@@..  ..::00@@@@oo    ::@@@@
@@oo..    ..oo@@@@@@@@00oo    ::@@@@00..      ..00@@@@oo    00@@
@@@@00::      ..........    ..00@@00::..::::..  ..@@@@@@::  ::@@
@@@@@@@@00::            ....oo@@@@::  ..00@@oo    00@@@@::  ..oo
..oo@@@@@@@@00oo::::oo0000@@@@@@@@..  ..@@@@@@::  ..oooo..      
    ..oo00@@@@@@@@@@@@@@@@@@@@@@@@::  ..00@@@@::      ........  
::      ..::00@@@@@@@@00oo::::00@@00..  ::00oo..    ..::oo0000oo
@@00::..      ..::......      ::@@@@00..      ..::00@@@@@@@@@@@@
@@@@@@00::..          ....    ..00@@@@..    ..oo@@@@@@@@@@0000@@
oo00@@@@@@@@oo::....::00@@oo  ..00@@@@..  ::00@@@@@@00::..  ....
  ..::00@@@@@@@@@@@@@@@@@@oo  ..@@@@oo    00@@@@00::
..      ::00@@@@@@@@@@@@oo..  oo@@@@..  ::@@@@00..    ..::oooo::
00oo..      ::oooooooo..    ::@@@@oo    oo@@@@..    ::00@@@@@@@@
@@@@@@oo..                ..00@@@@..  ::@@@@oo    oo@@@@@@@@@@@@
00@@@@@@00::........  ..::00@@@@::    00@@@@..  ::@@@@@@oo....::
  ..oo@@@@@@@@@@@@0000@@@@@@@@::    oo@@@@::    00@@@@::        
      ::@@@@@@@@@@@@@@@@@@@@::    ::@@@@00    ::@@@@oo    ..::..
oo..    ::oooooooooo000000::    ::@@@@@@::  ..@@@@00..  ::@@@@@@
@@@@oo                ....    ::@@@@@@::    oo@@@@::  ..00@@@@@@
@@@@@@::    ......          ..00@@@@::    ::@@@@oo    oo@@@@oooo
oo@@@@::  ::00@@@@@@00oo..  ::@@@@@@..    00@@@@..  ..@@@@oo
00@@00....00@@@@@@@@@@@@oo    oo@@@@::..::@@@@oo    oo@@@@..
@@@@::  ..@@@@oo::oo@@@@@@..  ..00@@@@00@@@@00..  ::@@@@00    ::
```


## Explanation

The Swift-Hohenberg Equation is given as  
```math
u'(t) = \varepsilon u - (\lambda^2 + \nabla^2)^2 \cdot u - u^3.
```
The right-hand side can be rewritten as a sum of a linear operator and a nonlinear operator on $u$:  
```math
u'(t) =  \quad   [\varepsilon - (\lambda^2 + \nabla^2)^2] \cdot u   \quad  +   \quad  -u^3  \quad = L(u) + N(u).
```
Starting with an initial condition $u(0)$ given, we can integrate with Euler's method:  
```math
u(t+1) = u(t) + dt * u'(t).
```
We make the method semi-implicit by representing $u'(t) = L(u(t+1)) + N(u(t))$, such that:  
```math
u(t+1) = u(t) + dt \cdot L(u(t+1))  +  dt \cdot N(u(t)).
```
Rearranging, we get the $u(t+1)$ terms on one side and the $u(t)$ terms on the other:  
```math
u(t+1) - dt \cdot L(u(t+1)) = u(t) +  dt \cdot N(u(t)).
```
Now, we take the Fourier Transform of both sides. The $\nabla^2\cdot$ operator becomes multiplication by $-k^2$ in the Fourier domain (for the 2D problem, this is a grid of frequency values). The resulting equation after applying an FFT is thus:  
```math
(1 - dt\cdot (\varepsilon - (\lambda^2 - k^2)^2)) * FFT[u(t+1)] = FFT[u(t) + dt\cdot N(u(t))].
```
For speed, we can precompute the linear operator array $Q = 1 - dt \cdot (\varepsilon - (\lambda^2 - k^2)^2)$ across the meshgrid of frequencies $k$. This gives the following in Fourier space:  
```math
FFT[u(t+1)] = FFT[u(t) + dt \cdot N(u(t))] / Q.
```
Thus, a full semi-implicit Euler time step update can be written as:  
```math
u(t+1) = iFFT\{FFT[u(t) + dt\cdot N(u(t))] / Q\}
```
In practice, this means:
1. Set initial grid `u(0)` as Gaussian noise
2. Precompute 2D meshgrid of `k` frequencies and the denominator term `Q = 1 - dt*(epsilon - (wavenum^2 - k^2)^2)`
3. Starting with current grid `u`, evaluate `u + dt * N(u)`
4. Compute the 2D FFT
5. Divide elementwise by the precomputed denominator term `Q`
6. Compute the 2D inverse FFT
7. Repeat steps 3-6 for desired number of time steps