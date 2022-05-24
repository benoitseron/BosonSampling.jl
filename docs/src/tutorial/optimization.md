# Optimization

An interesting question in BosonSampling is to find interferometers that of maximize certain properties.

<!-- For instance, consider a partition of all odd modes 1,3,... When using a `Fourier` interferometer, a strange phenomenon appears: -->

We provide the function `minimize_over_unitary_matrices()` which operates a conjugate gradient algorithm for the optimization over unitary matrices. It is implemented from [Conjugate gradient algorithm for optimization under unitary matrix constraint](https://doi.org/10.1016/j.sigpro.2009.03.015) from Traian Abrudan, Jan Eriksson, Visa Koivunen.
