```julia
# Define an input of 2 photons among 4 modes
i = Input{Bosonic}(first_modes(2,4))

# Define the interferometer
interf = RandHaar(4)

# Set the output measurement
o = FockSample()

# Create the event
ev = Event(i, o, interf)

# Simulate
sample!(ev)

# Scattershot boson sampling
scattershot_sampling(2,4)
```