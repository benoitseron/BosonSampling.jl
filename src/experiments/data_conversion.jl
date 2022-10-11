# convert CSV data into a vector of events

# conventions:
# the CSV file needs to have the same number of columns - make different files if using a different m

using DelimitedFiles
data = readdlm("data/loop_examples/indistinguishable/3-fold coincidences (2).csv", ',')

s = data[1,3]

replace(s, "(" => "")
replace(s, ")" => "")
replace(s, "\t" => "")



@withkw struct ExperimentalData
    data::Matrix{Int}
    n::Int
    m::Int
    n_samples::Union{Int, Nothing}
end



for line in 1:size(data)[1]

"""

"""
get_m(data::Matrix) = size(data)[2]
get_n_max(data::Matrix) = begin
    n_max = 0
    for line in 1:size(data)[1]
        n_max = maximum(n_max, sum(data[line]))
    end
    n_max
end

m = get_m(data)

n_subsets = 2
part = equilibrated_partition(m, n_subsets)
