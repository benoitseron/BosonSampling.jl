using JLD

cd("docs/publication/partitions/data/save")
data = load("number_samples_x.jld")

print(data["n_samples_array"])
