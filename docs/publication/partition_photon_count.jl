n = 10
m = 17
n_subsets = 4
input_state = first_modes(n,m)

interf = RandHaar(m)
i = Input{Bosonic}(input_state)
part = equilibrated_partition(m,n_subsets)
o = PartitionCountsAll(part)

ev = Event(i,o,interf)

compute_probability!(ev)

p_array = ev.proba_params.probability

# find the index of the most probable partition
max_proba_index = argmax(p_array.proba)
# display the corresponding photon counts
p_array.counts[max_proba_index]

# output:

# 3 in subset = [1, 2, 3, 4, 5]
# 2 in subset = [6, 7, 8, 9]
# 2 in subset = [10, 11, 12, 13]
# 3 in subset = [14, 15, 16, 17]
