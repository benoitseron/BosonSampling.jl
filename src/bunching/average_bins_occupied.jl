include("packages_loop.jl")

n_samples = 1000

begin
    n = 5
    m = n
    interf = RandHaar(m)
    T = OneParameterInterpolation
    x = 1.
    o = FockSample()
end

sample_number_modes_occupied(params::SamplingParameters) = number_modes_occupied(BosonSampling.sample!(params))

mean_number_modes_occupied(params::SamplingParameters, n_samples::Int) =  mean([sample_number_modes_occupied(params) for i in 1:n_samples])

mean_number_modes_occupied(params, n_samples)

x_array = 0:0.1:1
n_occ = []

for x in x_array
    params = SamplingParameters(n=n,m=m,interf = interf, T=T,x=x, o=o)
    set_parameters!(params)
    push!(n_occ,mean_number_modes_occupied(params, n_samples))
end

plot(x_array, n_occ)
