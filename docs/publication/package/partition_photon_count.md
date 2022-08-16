```julia
n = 10
m = 20

# experiment parameters
input_state = first_modes(n,m)
interf = RandHaar(m)
i = Input{Bosonic}(input_state)

# subset selection
s = Subset(first_modes(Int(m/2),m))
part = Partition(s)

# want to find all photon counting probabilities
o = PartitionCountsAll(part)

# define the event and compute probabilities
ev = Event(i,o,interf)
compute_probability!(ev)

# ... add some output
```