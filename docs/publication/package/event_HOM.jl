# Define the input with one particle in
# each input mode
r_i = ModeOccupation([1,1])
i = Input{Bosonic}(r_i)

# Define the unbalanced interferometer
B = BeamSplitter(1/sqrt(2))

# Define the output
r_f = ModeOccupation([1,1])
o = FockDetection(r_f)

# Create the event
ev = Event(i,o,B)
# Compute its probability to occur
compute_probability!(ev)
