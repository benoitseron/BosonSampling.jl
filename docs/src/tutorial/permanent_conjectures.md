# Permanent conjectures

Permanent and generalized matrix functions conjectures are linked to interesting and practical properties of boson samplers, as emphasized for instance by V. S. Shchesnovich in [Universality of Generalized Bunching and Efficient Assessment of Boson Sampling](https://arxiv.org/abs/1509.01561) as well as in the author's work [Boson bunching is not maximized by indistinguishable particles](https://arxiv.org/abs/2203.01306).

To search for new counter examples of a conjecture, one can implement a user-defined `search_function()`. For instance, `random_search_counter_example_bapat_sunder` searches for counter examples of the Bapat-Sunder conjecture (see also `violates_bapat_sunder`) in a brute-force manner, trying a different random set of matrices at each call.
One can then use

    search_until_user_stop(search_function)

which will iterate the function until you press `Ctrl+C` to interrupt the computation.

Another important conjecture is the permaent-on-top conjecture, disproved by V. S. Shchesnovich in [The permanent-on-top conjecture is false](https://www.sciencedirect.com/science/article/pii/S0024379515006631).
Special matrices related to this conjecture are given in this package such as the `schur_matrix(H)`, the general partial distinguishability function J(σ) implemented as `J_array`. From a matrix J, one can recover the density matrix of the internal states with `density_matrix_from_J`.
