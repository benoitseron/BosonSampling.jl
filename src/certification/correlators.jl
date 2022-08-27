# n = 8
# m = 8
# interf = RandHaar(m)
# TIn = Bosonic
# input_state = Input{TIn}(first_modes(n,m))

function C(i,j, interf::Interferometer, input_state::Input{T}) where {T <: Union{Bosonic, Distinguishable}}

    U = interf.U
    q = input_state.r.state

    coeff_d(i,j) = - sum([U[q[k], i] * U[q[k], j] * conj(U[q[k], i]) * conj(U[q[k], j])  for k in 1:n])

    correction_coeff_b(i,j) =  sum([k !=l ? (U[q[k], i] * U[q[l], j] * conj(U[q[l], i]) * conj(U[q[k], j])) : 0 for k in 1:n for l in 1:n])

    TIn = get_parametric_type(input_state)[1]

    if TIn == Bosonic
        return coeff_d(i,j) + correction_coeff_b(i,j)
    elseif TIn == Distinguishable
        return coeff_d(i,j)
    else
        error("not implemented")
    end

end

function C_dataset(interf::Interferometer, input_state::Input{T}) where {T <: Union{Bosonic, Distinguishable}}

    U = interf.U
    m = interf.m

    [C(i,j, interf, input_state) for i in 1:m for j in 1:m]

end

function correlators_nm_cv_s(interf::Interferometer, input_state::Input{T}) where {T <: Union{Bosonic, Distinguishable}}

    C_data = C_dataset(interf, input_state)

    moment(k) = mean(C_data.^k)

    moments = [moment(k) for k in 1:3]

    m = interf.m
    n = input_state.n

    nm = moments[1] * m^2 /n
    cv = sqrt( moments[2] -  moments[1]^2  ) /  moments[1]
    s = ( moments[3] - 3*  moments[1] * moments[2] -2 * moments[1]^3 )/(moments[2] -  moments[1]^2 )^(3/2)

    @warn "there seems to be a sign mistake!"
    (nm,cv,s)

end

n = 8
m = 8
interf = RandHaar(m)
TIn = Bosonic

events = generate_experimental_data(n_events = 100, n = n,m = m, interf = RandHaar(m), TIn = Bosonic)

function C(i,j, events::Vector{Event})

    count_in(i, ev::Event) = ev.output_measurement.s.state[i]

    counts_i = [count_in(i, ev) for ev in events]
    counts_j = [count_in(j, ev) for ev in events]

    mean(counts_i .* counts_j) - mean(counts_i) * mean(counts_j)

end
