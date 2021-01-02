# WARNING: This library's development has been suspended, possibly forever.
### Now the code is being rewritten in Rust programming language, since it's superior to C++ in many ways.
### All the future development will be conducted here: https://github.com/Graph-Donte-Crypto/AstroGraphicRust

# astro_cpp
Astrodynamics library in C++

Implemented features:
  + **Elliptic Kepler's equation solver**. Uses good initial guess and 3rd order Householder's method (about 40% faster than naive Newton's method).
  + **Lambert's problem solver** with multiple revolutions support. Based on Izzo's solver from PyKEP project (actually, this implementation is about 20% faster)

Future plans:
  + Gravity assists propagator (powered & unpowered)
  + r1, v1 --> r2, v2 (without intermediate conversion to orbital elements)
  + Optimizing trajectories involving Multiple Gravity Assists & Deep Space Maneuvers (MGADSM). It will employ Differential Evolution technique
  + Implementing GPU (CUDA, OpenCL, etc.) version of the optimizer
