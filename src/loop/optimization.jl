include("packages_loop.jl")

n = 2
sparsity = 3
m = sparsity * n
n_subsets = 3
n_subsets > 3 && m > 12 ? (@warn "may be slow") : nothing

x = 0.9


equilibrated_input(sparsity, m)

# photons are homogeneously distributed

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)

η_loss_lines = 1. * ones(m)
η_loss_bs = 1. * ones(m-1)



# η_0 = rand(m-1)
η_0 = 1/sqrt(2) * ones(m-1)
#η_0 = η_thermalization(m)


lower = zeros(length(η_0))
upper = ones(length(η_0))

#optimize(tvd_reflectivities, lower, upper, η_0)

sol = optimize(tvd_reflectivities, η_0, Optim.Options(time_limit = 60.0))

@show sol.minimizer

begin
    plot(sol.minimizer, label = "optim : $(-tvd_reflectivities(sol.minimizer))")
    plot!(η_thermalization(m), label = "thermalization : $(-tvd_reflectivities(η_thermalization(m)))")
    plot!(η_0, label = "initial : $(-tvd_reflectivities(η_0))")
    ylims!((0,1.3))
    ylabel!("η_i")
    xlabel!("loop pass")
end


### identical ###


tvd_one_reflectivity(η) = tvd_reflectivities(η[1] * ones(m-1))
η_0 = [0.5]

sol = optimize(tvd_one_reflectivity, η_0, Optim.Options(time_limit = 15.0))

@show sol.minimizer


### what happens with the thermalization pattern ###


(-tvd_reflectivities(η_thermalization(m)))

part = equilibrated_partition(m, n_subsets)
