include("packages_loop.jl")

# convert CSV data into a vector of events

# conventions:
# the CSV file needs to have the same number of columns - make different files if using a different m

using DelimitedFiles
using Dates

### working with csv ###

data = readdlm("data/loop_examples/indistinguishable/6-folds.csv", ',')

# needs to be in the form:
# 49, 1, 4, 5, 6, 7, 8, 10
# for finding 49 occurences of the mode list (1, 4, 5, 6, 7, 8, 10)

# conversion to be coded

# data
#
# s = data[1,:]
#
# replace(s, "(" => "")
# replace(s, ")" => "")
# replace(s, "\t" => "")
#


### run info ###

n = 4
sparsity = 2
m = sparsity * n

# x = 0.9

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)
η_loss_lines = nothing # 0.9 * ones(m)
η_loss_bs = nothing # 1. * ones(m-1)

η = 0.5 * ones(m-1)

params = LoopSamplingParameters(n=n, m=m, η = η, η_loss_bs = η_loss_bs, η_loss_lines = η_loss_lines, ϕ = ϕ)

### samples ###

# this is an example with fake samples, you will need to convert yours in the right format using a ModeList

samples = Vector{ThresholdModeOccupation}()
n_samples = 10

for i in 1:n_samples

    push!(samples, ThresholdModeOccupation(random_mode_list_collisionless(n,m)))

end

samples

### extra info on the run ###

extra_info = "this experiment was realised on... we faced various problems..."

### compiling everything in a single type structure ###

this_experiment = OneLoopData(params = params, samples = samples, extra_info = extra_info, name = "example")

### saving as a Julia format ###

save(this_experiment)

d = load("data/one_loop/example.jld")

loaded_experiment = d["data"]

loaded_experiment.samples
loaded_experiment.extra_info
