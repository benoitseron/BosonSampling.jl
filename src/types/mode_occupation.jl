"""
    ModeOccupation(state)

A list of the size of the number of modes `m`, with entry `j` of `state` being the number of photons in mode `j`. See also [`ModeList`](@ref).

    fields:
         - n::Int
         - m::Int
         - state::Vector{Int}
"""
@auto_hash_equals mutable struct ModeOccupation
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
    to_threshold!(mo::ModeOccupation)

Converts a `ModeOccupation` into threshold detection. Converts it into a `ThresholdModeOccupation` if not using the inplace version.
"""
function to_threshold(v::Vector{Int})

    [(mode >= 1 ? 1 : 0) for mode in v]

end

function to_threshold(mo::ModeOccupation)

    ThresholdModeOccupation(to_threshold(mo.state))

end

function to_threshold!(mo::ModeOccupation)

    mo.state = to_threshold(mo.state)

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
    n_detected::Int

    function ThresholdModeOccupation(ml::ModeList)

        state = convert(ModeOccupation, ml).state

        if !all(state[:] .>= 0)
            error("negative mode state")
        elseif !all(state[:] .<= 1)
            error("state can be at most one")
        else
            new(ml.m, state, sum(state))
        end
    end

    function ThresholdModeOccupation(mo::ModeOccupation)

            state = mo.state
            if !all(state[:] .>= 0)
                error("negative mode state")
            elseif !all(state[:] .<= 1)
                error("state can be at most one")
            else
                new(mo.m, state, sum(state))
            end

    end

    ThresholdModeOccupation(v::Vector{Int}) = ThresholdModeOccupation(ModeOccupation(v))




end

Base.sum(mo::ModeOccupation) = sum(mo.state)
Base.sum(mo::ThresholdModeOccupation) = sum(mo.state)

"""

    possible_threshold_detections_lossless(state::Vector{Int})

Returns a list of all possible states that can be obtained from `state`, the state resulting from a threshold detection.

"""
function possible_threshold_detections_lossless(n::Int, state::Vector{Int})

    # get all indexes with a photon

    indexes = findall(x->x==1, state)

    # get number of photons detected

    n_detected = sum(state)

    if n_detected == 0

        @show n, state
        error("no compatible lossless detections for this state")
        return []

    end

    if n_detected > n

        @show n, state
        error("incoherence: n_detected > n")

    end

    # if no loss, nothing to do 
    if n_detected == n
        return [state]
    end

    # finding the position of possible colliding photons

    @show n, n_detected
    mode_configs_colliding_photons = all_mode_configurations(n - n_detected, n_detected, only_photon_number_conserving = true)

    possible_states = []

    # generate a list of each possible configuration of colliding photons
    # add this configuration where a photon is detected in `state` (as labelled by `indexes`)

    for mode_config in mode_configs_colliding_photons

        # @show mode_config
     
        new_state = copy(state)
    
            # @show new_state
            # @show indexes
    
        for i in 1:length(indexes)
    
            new_state[indexes[i]] += mode_config[i]
    
        end
    
            # @show new_state
    
            push!(possible_states, new_state)
    
    end

    unique(possible_states)

end

"""

    possible_threshold_detections(n::Int, state::Vector{Int}; lossy::Bool = false)

Returns a list of all possible states that can be obtained from `state`, the state resulting from a threshold detection. If `lossy` is true, then the number of photons lost is also taken into account. This means computing all possible repartition of the lost photons among the last half of the output modes.

"""
function possible_threshold_detections(n::Int, state::Vector{Int}; lossy::Bool = false)

    if !lossy

        return possible_threshold_detections_lossless(n, state)

    else

        # find the number of modes

        m = length(state)

        # verify that the number of modes is even using @argcheck

        @argcheck iseven(m) "The number of modes must be even"

        # find half the number of modes as an Integer

        m_half = div(m,2)

        possible_states_lossless = possible_threshold_detections_lossless(n, state) 

        # get number of photons in the first half of each possible state

        n_not_lost_array = [sum(state[1:div(length(state),2)]) for state in possible_states_lossless]

        # get the number of lost photons for each possible state from n_not_lost_array

        n_lost_array = [n - n_not_lost for n_not_lost in n_not_lost_array]

        possible_states_lossy = []

        for possible_state_lossless in possible_states_lossless

            # @show "starting state: $possible_state_lossless"

            # get the number of lost photons for this possible state by counting the number of photon in the last half of this state

            n_lost = n - sum(possible_state_lossless[1:div(length(possible_state_lossless),2)])
            
            # generate possible mode configurations in the last half of the state for n_lost photons

            if n_lost == 0

                # @show "no lost photons"

                push!(possible_states_lossy, possible_state_lossless)

            else
                    
                mode_configs_lost_photons = all_mode_configurations(n_lost, m_half, only_photon_number_conserving = true)

                # add each possible configuration of lost photons to the possible state
                # push to possible_states_lossy

                # @show mode_configs_lost_photons

                for mode_config in mode_configs_lost_photons

                    new_state = copy(possible_state_lossless)

                    new_state[div(length(possible_state_lossless),2)+1:end] = mode_config

                    push!(possible_states_lossy, new_state)

                end

            end

        end
            
        return unique(possible_states_lossy)

    end

end


# write possible_treshold_detections for a ThresholdModeOccupation
# output a list of ModeOccupation

function possible_threshold_detections(n::Int,state::ThresholdModeOccupation; lossy = false)

    possible_states = possible_threshold_detections(n, state.state, lossy = lossy)

    possible_mode_occupations = []

    for state in possible_states

        push!(possible_mode_occupations, ModeOccupation(state))

    end

    possible_mode_occupations

end


# remove the lossy part of a ThresholdModeOccupation by removing the last half of the output modes and dividing m by two

function remove_lossy_part!(tmo::ThresholdModeOccupation)

    # check that the number of modes is even
    if tmo.m % 2 != 0
        error("The number of modes is not even")
    end

    tmo.m = tmo.m ÷ 2
    tmo.state = tmo.state[1:tmo.m]
    tmo.n_detected = sum(tmo.state)
    tmo
end

function remove_lossy_part(tmo::ThresholdModeOccupation)

    tmo_ = deepcopy(tmo)

    remove_lossy_part!(tmo_)

    tmo_

end

# write the same function for a ModeOccupation

function remove_lossy_part!(mo::ModeOccupation)

    # check that the number of modes is even
    if mo.m % 2 != 0
        error("The number of modes is not even")
    end

    mo.m = mo.m ÷ 2
    mo.state = mo.state[1:mo.m]
    mo
end

function remove_lossy_part(tmo::ModeOccupation)

    tmo_ = deepcopy(tmo)

    remove_lossy_part!(tmo_)

    tmo_

end


