@testset "sources" begin

    @test possible_inputs_loss([1,1,0,1], 3) == [[1, 1, 0, 1], [0, 1, 0, 1], [1, 0, 0, 1], [1, 1, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0]]

    # test that possible_inputs_loss([1,0,1], 3) throws an ArgumentError

    @test_throws ArgumentError possible_inputs_loss([1,0,1], 3)


    @test input_probability([1,1,0,1], [1,1,0,1], QuantumDot(efficiency=1)) == 1
    @test input_probability([1,1,0,1], [1,0,0,1], QuantumDot(efficiency=1)) == 0
    @test input_probability([1,1,0,1], [1,1,0,1], QuantumDot(efficiency=0)) == 0

    # now for imperfect source probabilities


    n = 3
    interf = Fourier(n)
    o = FockDetection(ModeOccupation([1,1,0]))

    source = QuantumDot(efficiency = 0.4)

    params = SamplingParameters(n = n, interf = interf, o = o)

    set_parameters!(params)

    # test that compute_probability!(params) fails

    @test_throws ErrorException compute_probability!(params)
    # this doesn't work as there are some lost photons in this lossless interferometer

    compute_probability_imperfect_source(params, source) = 0.09600000000000009
end