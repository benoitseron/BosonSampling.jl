# Samplers

This tutorial gives some examples of the usage for the samplers. That is, from an
[`Input`](@ref) configuration and an [`Interferometer`](@ref) we simulate the
resulting configuration.

## Bosonic sampler

We present here the general syntax to simulate `n=4` indistinguishable photons among
`m=16` modes. To do so, we first need to define our [`Bosonic`](@ref) input with
randomly placed photons

    julia> n = 4;

    julia> m = n^2;

    julia> my_input = Input{Bosonic}(ModeOccupation(random_occupancy(n,m)))
     Type:Input{Bosonic}
     r:state = [0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1]
     n:4
     m:16
     G:GramMatrix{Bosonic}(4, ComplexF64[1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im], nothing, nothing, OrthonormalBasis(nothing))
     distinguishability_param:nothing

and we use a random interferometer

    julia> my_interf = RandHaar(m)
    Interferometer :

     Type : RandHaar
     m : 16

and then call [`cliffords_sampler`](@ref) to run the simulation

    julia> res = cliffords_sampler(input=my_input, interf=my_interf)
    4-element Vector{Int64}:
     2
     8
     15
     16

The output vector of length `n` tells us which of the output modes contain a photon. One can have a schematic look at the input/output configurations:

    julia> visualize_sampling(my_input, res)


## Noisy sampler

A boson sampling experiment can also be simulated by adding noise either from the interferometer or directly by playing with partial distinguishability. As before, one creates an input of particles that are not completely indistinguishable from [`OneParameterInterpolation`](@ref)

    julia> my_distinguishability_param = 0.7;

    julia> my_mode_occupation = ModeOccupation(random_occupancy(n,m));

    julia> my_input = Input{OneParameterInterpolation}(my_mode_occupation, my_distinguishability_param);

and still using `my_interf` with some loss `η=0.7`, one simulates our noisy boson
sampling experiment  with

    julia> res = noisy_sampler(input=my_input, reflectivity=η, interf=my_interf)
    3-element Vector{Int64}:
     5
     7
     11

where we have lost one particle as `length(res)=3`, meaning that only three output modes
are populated by one photon.

## Classical sampler

Finally, we repeat the steps to simulate fully distinguishable particles by using
[`classical_sampler`](@ref)

    julia> my_mode_occupation = ModeOccupation(random_occupancy(n,m));

    julia> my_input = Input{Distinguishable}(my_mode_occupation);

    julia> my_interf = RandHaar(m);

    julia> res = classical_sampler(input=my_input, interf=my_interf)
     16-element Vector{Int64}:
     0
     0
     0
     1
     ⋮
     0
     1
     0
