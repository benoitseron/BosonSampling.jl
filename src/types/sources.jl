"""
    abstract type Source 

Holds types simulating realistic experimental sources.
"""
abstract type Source end

"""

    mutable struct QuantumDot <: Source

Type holding a model of a QuantumDot. The point is to simulate its non deterministic photon generation, to know what kind of input is sent through the interferometer.

This is held in through the field `efficiency`, the probability that a photon is generated if we ask one in this position. The probabilities are assumed to be IID.
"""
@with_kw mutable struct QuantumDot <: Source
    efficiency::Real = 1. # probability that a photon is generated if we ask one in this position
end



"""

    possible_inputs_loss(input_no_loss::Vector{Int}, lost::Int)

Generates all possible input patterns from an ideal `input_no_loss` and with up to `lost` photons

Example:

    possible_inputs_loss([1,1,0,1], 3)

    8-element Vector{Any}:
    [1, 1, 0, 1]
    [0, 1, 0, 1]
    [1, 0, 0, 1]
    [1, 1, 0, 0]
    [0, 0, 0, 1]
    [0, 1, 0, 0]
    [1, 0, 0, 0]
    [0, 0, 0, 0]
"""
function possible_inputs_loss(input_no_loss::Vector{Int}, lost::Int)

    @argcheck is_collisionless(input_no_loss) "input_no_loss has more than one photon per mode, not implemented"

    @argcheck lost <= sum(input_no_loss) "lost photons is larger than the number of photons in input_no_loss"

    m = length(input_no_loss)

    # get the indexes with a photon

    indexes = findall(x->x==1, input_no_loss)

    # generate a list of all possible indexes, deleting one at a time, up to `lost`

    photon_deleter(l_lost) = unique(collect(permutations([i > l_lost for i in 1:length(indexes)])))

    function indexes_remaining(l_lost)  
        
        indexes_deletion = [indexes .* photon_deleter(l_lost)[i] for i in 1:length(photon_deleter(l_lost))]

        [indexes_deletion[el][findall(x-> x != 0, indexes_deletion[el])]
    for el in 1:length(indexes_deletion)]

    end

    # convert remaining indexes to remaining photons 

    function possible_inputs(l_lost)

        result = []

        for indexes in indexes_remaining(l_lost)

            input = zeros(Int,m)

            for index in indexes

                input[index] = 1

            end

            push!(result, input)

        end

        unique(result)
    end


    all_inputs = []

    for l_lost in 0:lost

        append!(all_inputs, possible_inputs(l_lost))

    end

    all_inputs

end


"""

    input_probability(input_ideal::Vector{Int}, input_real::Vector{Int}, source::QuantumDot)


Computes the probability of obtaining `input_real` from `input_ideal` for a simple quantum dot source.
"""
function input_probability(input_ideal::Vector{Int}, input_real::Vector{Int}, source::QuantumDot)

    # compute number of modes

    m = length(input_ideal)

    # compute number of photons

    n = sum(input_ideal)

    # compute number of photons in the real input

    n_real = sum(input_real)

    # assert that the number of modes is the same

    @argcheck m == length(input_real)

    # compute the probability to generate n_real photons when n photons are generated, each one with a probability given by the quantum dot `efficiency` through a binomial distribution

    proba = binomial(n, n_real) * source.efficiency^n_real * (1-source.efficiency)^(n-n_real)

end

