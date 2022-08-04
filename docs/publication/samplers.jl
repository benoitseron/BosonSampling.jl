# Define an input of 3 photons among 3 modes
i = Input{Bosonic}(first_modes(3,3))

# Define the interferometer
interf = RandHaar(3) 

# Set the output measurement
o = FockSample()

# Create the event
ev = Event(i, o, interf)

# Simulate
sample!(ev)
