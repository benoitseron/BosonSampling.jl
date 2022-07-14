"""
    virtual_interferometer_uniform_loss(real_interf, η)

Simulates a simple, uniformly lossy interferometer: take a 2m*2m interferometer and introduce beam splitters in front with transmittance `η`. Returns the corresponding virtual interferometer.
"""
function virtual_interferometer_uniform_loss(real_interf, η)

    # define a (2m) * (2m) virtual interferometer V

    # connect the remaining branches of the input beam splitters
    # (of course this could be made different but this is the simplest,
    # no thought case)

    U = real_interf.U

    m = size(U,1)
    V = Matrix{eltype(U)}(I, (2m, 2m))

    # loss beam splitters
    for i in 1:m
        V *= beam_splitter_modes(in_up = i, in_down = i+m, out_up = i, out_down = m+i, transmission_amplitude = η, n = 2m)
    end

    U_large = Matrix{eltype(U)}(I, (2m, 2m))
    U_large[1:m,1:m] = U
    V *= U_large

    #V = V[1:2m, 1:2m] # disregard virtual mode to connect second input branch of beam splitters
    ############ this must clearly not be good for unitarity

    V = copy(transpose(V))

    UserDefinedInterferometer(V)

end


"""
    to_lossy(s::Subset)

Transforms a subset of size `m` into a `2m` one with the last `m` modes being a empty (the environment modes).
"""
function to_lossy(s::Subset)
    subset = s.subset
    new_subset = zeros(Int, 2*length(subset))
    new_subset[1:length(subset)] .= subset
    Subset(new_subset)
end

"""
    to_lossy(part::Partition)

Transforms a partition of size `m` into a `2m` one with the last `m` modes being a subset (the environment modes).
"""
function to_lossy(part::Partition)

    environment = Subset(last_modes(n,2m))

    new_subsets =  Vector{Subset}()

    for subset in part.subsets
        push!(new_subsets, to_lossy(subset))
    end

    push!(new_subsets, environment)

    Partition(new_subsets)

end
