# BosonSampling

This project implements standard and scattershot BosonSampling in Julia, including boson samplers and certification tools.

## Conventions

### Basic conventions :

  a_j -> \sum_j U_jk b_k

  so the rows correspond to the input, columns to the output, that is : the probability that a single goes from j to k is |U_jk|^2

  By default, we will use Tichy's conventions as in
  https://arxiv.org/abs/1312.4266
  Input vector = r or input_state
  Output vector = s or output_state (output_measurement as well)
  interferometer matrix = U
  interferometer matrix M of tichy (with rows corresponding to the input,...) = scattering_matrix
  dimension of the matrix = m (size of the interferometer) or (n in previous code
  or where the number of photons is irrelevant or called number_photons)
  number of modes occupied = n, number_photons

### Conventions regarding Julia:

  Unlike most languages, the counting goes from 1,2,3... instead of starting at
  zero as 0,1,2,...

  We write matrices as Matrix(Line, Column). Parts of the program were written
  with the opposite convention, thus there may be mistakes although most of them
  should be gone. (I believe Julia's convention changed in a version update).

### Warning for Gram matrices :

  typically it is defined as (<phi_i|phi_j>); i,j = 1:n the label of the photon
  if photons are swapped, you need to enter another distinguishability matrix
  swapped accordingly

### Warning about precision :

  in EventProbability :

  precision is set to machine precision eps() when doing non-randomised methods
  although it is of course larger and this should be implemented
  with permanent approximations, see for instance
  https://arxiv.org/abs/1904.06229

### Distances :

  Beware of the different TVD conventions (1/2 in front or not)
