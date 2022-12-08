"""
    ModeOccupation(state)

A list of the size of the number of modes `m`, with entry `j` of `state` being the number of photons in mode `j`. See also [`ModeList`](@ref).

    fields:
         - n::Int
         - m::Int
         - state::Vector{Int}
"""
@auto_hash_equals struct ModeOccupation
    n::Int
    m::Int
    state::Vector{Int}
    ModeOccupation(state) = all(state[:] .>= 0) ? new(sum(state), length(state), state) : error("negative photon counts")
end

Base.show(io::IO, i::ModeOccupation) = print(io, "state = ", i.state)

"""
        number_modes_occupied(mo::ModeOccupation)

Number of modes having at least a photon.
"""
number_modes_occupied(mo::ModeOccupation) = sum(to_threshold(mo).state)


"""

    :+(s1::ModeOccupation, s2::ModeOccupation)
    :+(s1::ModeOccupation, s2::Vector{Int})
    :+(s2::Vector{Int}, s1::ModeOccupation)

Adds two mode occupations, for instance
s1 = ModeOccupation([0,1])
s2 = ModeOccupation([1,0])

(s1+s2).state == [1,1]

Also works with just a vector and a mode occupation.
"""
Base.:+(s1::ModeOccupation, s2::ModeOccupation) = begin
    return ModeOccupation(s1.state + s2.state)
end

Base.:+(s1::ModeOccupation, s2::Vector{Int}) = begin

        @argcheck size(s1.state) == size(s2) "incompatible sizes"
    return ModeOccupation(s1.state + s2)
end

Base.:+(s2::Vector{Int}, s1::ModeOccupation) = begin
    return s1 + s2
end

"""
        Base.zeros(mo::ModeOccupation)

Returns a `ModeOccupation` similar to the input but with a state made of zeros.
"""
function Base.zeros(mo::ModeOccupation)

    physical_state = mo.state
    state = zeros(eltype(physical_state), size(physical_state))
    ModeOccupation(state)

end




"""
    to_threshold(v::Vector{Int})
    to_threshold(mo::ModeOccupation)

Converts a `ModeOccupation` into threshold detection.
"""
function to_threshold(v::Vector{Int})

    [(mode >= 1 ? 1 : 0) for mode in v]

end

function to_threshold(mo::ModeOccupation)

    ModeOccupation(to_threshold(mo.state))

end


"""
        Base.cat(s1::ModeOccupation, s2::ModeOccupation)

Concatenates two `ModeOccupation`.
"""
function Base.cat(s1::ModeOccupation, s2::ModeOccupation)

    ModeOccupation(vcat(s1.state, s2.state))

end



"""
    ModeList(state)
    ModeList(state,m)

Contrasting to [`ModeOccupation`](@ref) this list is of size `n`, the number of photons. Entry `j` is the index of the mode occupied by photon `j`.

This can also be used just to select modes for instance.

See also [`ModeOccupation`](@ref).

    fields:
        - n::Int
        - m::Union{Int, Nothing}
        - modes::Vector{Int}
"""
@auto_hash_equals struct ModeList
    n::Int
    m::Union{Int, Nothing}
    modes::Vector{Int}

    ModeList(modes::Vector{Int}) = ModeList(modes, nothing)

    # all(modes[:] .>= 1) ? new(length(modes), nothing, modes) : error("modes start at one")

        function ModeList(modes::Vector{Int}, m)

                if all(modes[:] .>= 1) && (m == nothing ? true : all(modes[:] .<= m))
                    new(length(modes), m, modes)
                else
                    error("incoherent or invalid mode inputs")
                end
        end

        ModeList(mode::Int, m = nothing) = ModeList([mode],m)

end

"""
        is_compatible(target_modes_in::ModeList, target_modes_out::ModeList)

Checks compatibility of `ModeList`s.
"""
function is_compatible(target_modes_in::ModeList, target_modes_out::ModeList)

        if target_modes_in == target_modes_out
                return true
        else
                @argcheck target_modes_in.n == target_modes_out.n
                @argcheck target_modes_in.m == target_modes_out.m
                true
        end

end

function Base.convert(::Type{ModeOccupation}, ml::ModeList)

        if ml.m == nothing
                error("need to give m")
        else
                state = zeros(Int, ml.m)

                for mode in ml.modes
                        state[mode] += 1
                end

                return ModeOccupation(state)
        end
end

function Base.convert(::Type{ModeList}, mo::ModeOccupation)

        mode_list = Vector{Int}()

        for (mode, n_in) in enumerate(mo.state)
            if n_in > 0
                for photon in 1:n_in
                    push!(mode_list, mode)
                end
            end
        end

        ModeList(mode_list, mo.m)

end

at_most_one_photon_per_bin(state) = all(state[:] .<= 1)
at_most_one_photon_per_bin(r::ModeOccupation) = at_most_one_photon_per_bin(r.state)

isa_subset(subset_modes::Vector{Int}) = (at_most_one_photon_per_bin(subset_modes) && sum(subset_modes) != 0)
isa_subset(subset_modes::ModeOccupation) = isa_subset(subset_modes.state)

"""
    first_modes(n::Int, m::Int)

Create a [`ModeOccupation`](@ref) with `n` photons in the first sites of `m` modes.
"""
first_modes(n::Int,m::Int) = n<=m ? ModeOccupation([i <= n ? 1 : 0 for i in 1:m]) : error("n>m")

first_modes_array(n::Int,m::Int) = first_modes(n,m).state

"""
    last_modes(n::Int, m::Int)

Create a [`ModeOccupation`](@ref) with `n` photons in the last sites of `m` modes.
"""
last_modes(n::Int,m::Int) = n<=m ? ModeOccupation([i > m-n ? 1 : 0 for i in 1:m]) : error("n>m")

last_modes_array(n::Int,m::Int) = last_modes(n,m).state

equilibrated_input(sparsity, m) = ModeOccupation([((i-1) % sparsity) == 0 ? 1 : 0 for i in 1:m])

"""
    mutable struct ThresholdModeOccupation

Holds threshold detector state. Example

    ThresholdModeOccupation(ModeList([1,2,4], 4))

"""
@auto_hash_equals mutable struct ThresholdModeOccupation

    m::Int
    state::Vector{Int}

    function ThresholdModeOccupation(ml::ModeList)

        state = convert(ModeOccupation, ml).state

        if !all(state[:] .>= 0)
            error("negative mode state")
        elseif !all(state[:] .<= 1)
            error("state can be at most one")
        else
            new(ml.m, state)
        end
    end

    function ThresholdModeOccupation(mo::ModeOccupation)

            state = mo.state
            if !all(state[:] .>= 0)
                error("negative mode state")
            elseif !all(state[:] .<= 1)
                error("state can be at most one")
            else
                new(mo.m, state)
            end

    end

    ThresholdModeOccupation(v::Vector{Int}) = ThresholdModeOccupation(ModeOccupation(v))




end

Base.sum(mo::ModeOccupation) = sum(mo.state)
Base.sum(mo::ThresholdModeOccupation) = sum(mo.state)

