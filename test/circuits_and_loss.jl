@testset "circuits and loss" begin

### 2d HOM with loss example ###

n = 1
m = 1

function lossy_line_example(η_loss)

    circuit = LossyCircuit(1)
    interf = LossyLine(η_loss)
    target_modes = [1]

    add_element_lossy!(circuit, interf, target_modes_in = target_modes)
    circuit

end

lossy_line_example(0.9)

transmission_amplitude_loss_array = 0:0.1:1
output_proba = []

i = Input{Bosonic}(to_lossy(first_modes(n,m)))
o = FockDetection(to_lossy(first_modes(n,m)))

for transmission in transmission_amplitude_loss_array

    ev = Event(i,o, lossy_line_example(transmission))
    compute_probability!(ev)
    push!(output_proba, ev.proba_params.probability)
end

@test output_proba ≈ [0.0, 0.010000000000000002, 0.04000000000000001, 0.09, 0.16000000000000003, 0.25, 0.36, 0.48999999999999994, 0.6400000000000001, 0.81, 1.0]

### 2d HOM with loss example ###

n = 2
m = 2
i = Input{Bosonic}(first_modes(n,m))
o = FockDetection(ModeOccupation([2,0])) # detecting bunching, should be 0.5 in probability if there was no loss
transmission_amplitude_loss_array = 0:0.1:1
output_proba = []

function lossy_bs_example(η_loss)

    circuit = LossyCircuit(2)
    interf = LossyBeamSplitter(1/sqrt(2), η_loss)
    target_modes = [1,2]

    add_element_lossy!(circuit, interf, target_modes_in = target_modes)
    circuit

end

for transmission in transmission_amplitude_loss_array

    ev = Event(i,o, lossy_bs_example(transmission))
    compute_probability!(ev)
    push!(output_proba, ev.proba_params.probability)
end

@test output_proba ≈ [0.0, 5.0000000000000016e-5, 0.0008000000000000003, 0.004049999999999998, 0.012800000000000004, 0.031249999999999993, 0.06479999999999997, 0.12004999999999996, 0.20480000000000007, 0.32805, 0.4999999999999999]

### building the loop ###

n = 3
m = n
η = 1/sqrt(2) .* ones(m-1)
# 1/sqrt(2) .* [1,0] #ones(m-1) # see selection of target_modes = [i, i+1] for m-1
# [1/sqrt(2), 1] #1/sqrt(2) .* ones(m-1) # see selection of target_modes = [i, i+1] for m-1

# η_loss = 1. .* ones(m-1)

circuit = LosslessCircuit(m) #LossyCircuit(m)

for mode in 1:m-1

    interf = BeamSplitter(η[mode]) #LossyBeamSplitter(reflectivities[mode], η_loss[mode])
    target_modes_in = [mode, mode+1]
    target_modes_out = [mode, mode+1]
    add_element!(circuit, interf, target_modes_in = target_modes_in, target_modes_out = target_modes_out)

end



i = Input{Bosonic}(first_modes(n,m))

#outputs compatible with two photons top mode
o1 = FockDetection(ModeOccupation([2,1,0]))
o2 = FockDetection(ModeOccupation([2,0,1]))

o_array = [o1,o2]

p_two_photon_first_mode = 0

for o in o_array
    ev = Event(i,o, circuit)
    compute_probability!(ev)
    p_two_photon_first_mode += ev.proba_params.probability
end

@test p_two_photon_first_mode ≈ 0.5

o3 = FockDetection(ModeOccupation([3,0,0]))
ev = Event(i,o3, circuit)
compute_probability!(ev)

@test ev.proba_params.probability ≈ 0.






end