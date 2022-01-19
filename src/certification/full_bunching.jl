# implements the certification algorithm of Valery S. Shchesnovich
# from https://arxiv.org/abs/1509.01561

function choose_subset_size(;m,n, minimal_bunching_proba = 0.25)

    """finds the optimal number of modes in the subset by maximizing the
    haar-averaged bunching probability,

    returns (ideal subset size, average full bunching probability ratio bosonic/dist)"""

    max_ratio = 0
    subset_size_max_ratio = nothing

    for subset_size in 1:m-1

        proba_dist, proba_bosonic = subset_expectation_value(subset_size,n,n,m)

        if proba_bosonic >= minimal_bunching_proba && !(proba_dist ≈ 0)  && proba_bosonic/proba_dist > max_ratio
            max_ratio = proba_bosonic/proba_dist
            subset_size_max_ratio = subset_size
        end
    end

    subset_size_max_ratio, max_ratio

end

choose_subset_size(m = 20,n = 4)
