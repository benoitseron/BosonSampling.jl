# Conventions

Photon creation operators are changed as
``a_j -> \sum_j U_jk b_k``
when going through the interferometer ``\op{U}``.
Thus, the rows correspond to the input, columns to the output, that is: the probability that a single goes from `j` to `k` is ``|U_jk|^2``

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

See for Tichy: https://arxiv.org/abs/1410.7687, and for Shchesnovich: https://arxiv.org/abs/1410.1506

## Operators change

* Tichy: ```\hat{A}_{j,|\Phi\rangle}^{\dagger} \rightarrow \hat{U} \hat{A}_{j,|\Phi\rangle}^{\dagger} \hat{U}^{-1}=\sum_{k=1}^{m} U_{j, k} \hat{B}_{k,|\Phi\rangle}^{\dagger}```

* Shchesnovich: ```A spatial unitary network can be defined by an unitary transformation between input $a_{k, s}^{\dagger}(\omega)$ and output $b_{k, s}^{\dagger}(\omega)$ photon creation operators, we set $a_{k, s}^{\dagger}(\omega)=\sum_{l=1}^{M} U_{k l} b_{l, s}^{\dagger}(\omega)$, where $U_{k l}$ is the unitary matrix describing such an optical network.```

This means that if Tichy uses U, then Shchesnovich has in place U^†.

## Scattering matrix

* Tichy: ```Using the mode assignment list $\vec{d}(\vec{s})=\left(d_{1}, \ldots d_{n}\right)$ [49], which indicates the mode in which the $j$ th particle resides, the effective scattering matrix becomes
$$
M=U_{\vec{d}(\vec{r}), \vec{d}(\vec{s})}
$$
where our convention identifies the $j$ th row (column) with the $j$ th input (output) mode, as illustrated in Fig. 1)(a).```

* Shchesnovich: identical

## Bunching

The H-matrix follows a convention different from that of Valery Shchesnovich: `H_{a,b} = \sum _{l \in \mathcal{K}} U_{l,a} U_{l,b}^{*}`, see [H_matrix](@ref).


## Conventions regarding Julia:

  Unlike most languages, the counting goes from 1,2,3... instead of starting at
  zero as 0,1,2,...

## Gram matrices :

  Gram matrices are defined as ``(<\phi_i|\phi_j>); i,j = 1:n``. This means that if the label of the photons are swapped, you need to enter another distinguishability matrix with swapped labels accordingly. See [GramMatrix](@ref).

## Warning about precision :

  in [EventProbability](@ref) :

  precision is set to machine precision `eps()` when doing non-randomised methods
  although it is of course larger and this should be implemented
  with permanent approximations, see for instance
  https://arxiv.org/abs/1904.06229

## Distances :

  Beware of the different TVD conventions (1/2 in front or not). See [tvd](@ref) for instance.
