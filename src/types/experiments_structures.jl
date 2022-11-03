@auto_hash_equals mutable struct ThresholdModeOccupation

    m::Int
    clicks::Vector{Int}

    function ThresholdModeOccupation(ml::ModeList)

        clicks = convert(ModeOccupation, ml).state

        if !all(clicks[:] .>= 0)
            error("negative mode clicks")
        elseif !all(clicks[:] .<= 1)
            error("clicks can be at most one")
        else
            new(ml.m, clicks)
        end
    end

end

# example ThresholdModeOccupation(ModeList([1,2,4], 4))

@with_kw mutable struct OneLoopData

    params::LoopSamplingParameters
    samples::Vector{ThresholdModeOccupation}
    extra_info

end

extra_info = "this experiment was realised on... we faced various problems..."

typeof(extra_info)

OneLoopData()


n = 10
sparsity = 2
m = sparsity * n

# x = 0.9

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)
η_loss_lines = nothing # 0.9 * ones(m)
η_loss_bs = nothing # 1. * ones(m-1)

η = 0.5 * ones(m-1)

params = LoopSamplingParameters(n=n, m=m, η = η, η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, ϕ = ϕ)

samples = Vector{ThresholdModeOccupation}()
n_samples = 10

for i in 1:n_samples

    push!(samples, ThresholdModeOccupation(random_mode_list_collisionless(n,m)))

end

samples
