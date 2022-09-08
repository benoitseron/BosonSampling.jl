```julia
n = 20
m = 400

# Define an input of 20 photons among 400 modes
i = Input{Bosonic}(first_modes(n,m))

# Define the interferometer
interf = RandHaar(m)

# Set the output measurement
o = FockSample()

# Create the event
ev = Event(i, o, interf)

# Simulate
sample!(ev)
# output:
# state = [0,1,0,...]
```