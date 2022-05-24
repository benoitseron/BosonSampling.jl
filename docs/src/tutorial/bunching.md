# Bunching

Boson bunching is at the heart of many quantum phenomena, and this package has multiple functions related to it in the context of linear optics.

Given an interferometer `interf`, the probability to find all photons of a given input `i` (with a general state of distinguishability) in a subset `subset` of the output modes is given by

    full_bunching_probability(interf::Interferometer, i::Input, subset_modes::Subset)

At the heart of this forumla lies the `H_matrix(interf::Interferometer, i::Input, subset_modes::ModeOccupation)`, describing the bunching properties of an interferometer and subset (see [Boson bunching is not
maximized by indistinguishable particles](https://arxiv.org/abs/2203.01306)).

Although inefficient, we also provide a check function to evaluate by direct summation the bunching probabilities for `Bosonic` inputs

    bunching_probability_brute_force_bosonic(interf::Interferometer, i::Input, subset_modes::ModeOccupation)
in order to check the implementation of the above functions.
