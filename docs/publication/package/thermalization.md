```julia
using BosonSampling
using Plots; pyplot()

n = 20
trunc = 10
T = OneParameterInterpolation
for x in 0:0.01:1
    i =
    Input{T}(first_modes(n,n),x)
    set1 = [0 for i in 1:n]
    set1[1] = 1
    part = Partition([Subset(set1)])
    F = Fourier(n)
    (idx,pdf) =
    compute_probabilities_partition(F,part,i)
end
```