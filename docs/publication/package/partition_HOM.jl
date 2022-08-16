i = Input{Bosonic}(first_modes(2,2))

set1 = Subset([1,0])
set2 = Subset([0,1])
interf = Fourier(2)
part = Partition([set1,set2])

(idx,pdf) =
compute_probabilities_partition(interf,part,i)
