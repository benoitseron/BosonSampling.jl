```julia
n = 25
m = 400

# Experiment parameters
input_state = first_modes(n,m)
interf = RandHaar(m)
i = Input{Bosonic}(input_state)

# Subset selection
s = Subset(first_modes(Int(m/2),m))
part = Partition(s)

# Want to find all photon counting probabilities
o = PartitionCountsAll(part)

# Define the event and compute probabilities
ev = Event(i,o,interf)
compute_probability!(ev)
# About 30s execution time on a single core
#
# output:
#
# 0 in subset = [1, 2,..., 200]
# p = 4.650035467008141e-8
# ...
```