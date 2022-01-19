### scattering ###

m = 5
n = 3

interf = RandHaar(m)
i = Input{PartDist}(first_modes(n,m), GramMatrix{PartDist}(n, rand_gram_matrix(n)))
o = OutputMeasurement{FockDetection}(first_modes(n,m))
ev = Event(i,o,interf)

compute_probability!(ev)
