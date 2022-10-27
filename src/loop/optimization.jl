include("packages_loop.jl")

n = 10
sparsity = 2
m = sparsity * n
n_subsets = 2
n_subsets > 3 && m > 12 ? (@warn "may be slow") : nothing

equilibrated_input(sparsity, m) = ModeOccupation([((i-1) % sparsity) == 0 ? 1 : 0 for i in 1:m])

# photons are homogeneously distributed

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)

η_loss_lines =1. * ones(m)
η_loss_bs = 1. * ones(m-1)

function tvd_reflectivities(η)

    if all(between_one_and_zero.(η))

        params = LoopSamplingParameters(n=n, m=m, η = η, η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, ϕ = ϕ)

        psp_b = convert(PartitionSamplingParameters, params)
        psp_d = convert(PartitionSamplingParameters, params) # a way to act as copy

        psp_b.mode_occ = equilibrated_input(sparsity, m)
        psp_d.mode_occ = equilibrated_input(sparsity, m)

        part = equilibrated_partition(m, n_subsets)

        psp_b.part = part
        psp_d.part = part

        psp_b.T = OneParameterInterpolation
        psp_b.x = 0.9
        set_parameters!(psp_b)

        psp_d.T = Distinguishable
        set_parameters!(psp_d)

        compute_probability!(psp_b)
        compute_probability!(psp_d)

        pdf_bos = psp_b.ev.proba_params.probability.proba
        pdf_dist = psp_d.ev.proba_params.probability.proba

        @show tvd(pdf_bos,pdf_dist)
        #display(plot(η))
        return -tvd(pdf_bos,pdf_dist)
    else
        println("invalid reflectivity")
        return 100
    end

end


# η_0 = rand(m-1)
# η_0 = 1/sqrt(2) * ones(m-1)
η_0 = η_thermalization(m)


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
