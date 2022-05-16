# About

This project implements standard and scattershot BosonSampling in Julia, including boson samplers and certification tools.

## Conventions

### Basic conventions :
  Photon creation operators are changed as ``a_j \right_arrow \sum_j U_{jk} b_k``
  when going through the interferometer ``\op{U}``.
  Thus, the rows correspond to the input, columns to the output, that is: the probability that a single goes from `j` to `k` is ``|U_{jk}|^2``

  This is the conventions used by most people. Let us warn that Valery Shchesnovich uses a convention that is incompatible: ``\op{U}`` needs to be changed to ``\op{U}^\dagger``. (And likewise defines the Gram matrix as the transpose of ours, see below.)

  By default, we will use [Tichy's conventions](https://arxiv.org/abs/1312.4266)
  * Input vector = `r` or `input_state`
  * Output vector = `s` or `output_state`
  * For detection that if not just an event `output_measurement`
  * interferometer matrix = `U`
  * interferometer matrix `M` of Tichy (with rows corresponding to the input,...) = `scattering_matrix`
  * dimension of the interferometer = `m` (size of the interferometer) or (`n` in previous code
  or where the number of photons is irrelevant or called number_photons)
  * number of modes occupied = `n`, `number_photons`

### Bunching

The H-matrix follows a convention different from that of Valery Shchesnovich: ``H_{a,b} = \sum _{l \in \mathcal{K}} U_{l,a} U_{l,b}^{*}``.


### Conventions regarding Julia:

  Unlike most languages, the counting goes from 1,2,3... instead of starting at
  zero as 0,1,2,...

### Gram matrices :

  Gram matrices are defined as ``(<\phi_i|\phi_j>); i,j = 1:n``. This means that if the label of the photons are swapped, you need to enter another distinguishability matrix with
  swapped labels accordingly.

### Warning about precision :

  in EventProbability :

  precision is set to machine precision `eps()` when doing non-randomised methods
  although it is of course larger and this should be implemented
  with permanent approximations, see for instance
  https://arxiv.org/abs/1904.06229

### Distances :

  Beware of the different TVD conventions (1/2 in front or not)
