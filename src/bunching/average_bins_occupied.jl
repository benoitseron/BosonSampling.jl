include("packages_loop.jl")

n = 3
m = n
interf = Fourier(m)
T = OneParameterInterpolation
x = 1.
o = FockSample()

params = SamplingParameters(n=n,m=m,interf = interf, T=T,x=x, o=o)

set_parameters!(params)

BosonSampling.sample!(params)
