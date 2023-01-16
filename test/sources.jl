@testset "sources" begin

    @test possible_inputs_loss([1,1,0,1], 3) == [[1, 1, 0, 1], [0, 1, 0, 1], [1, 0, 0, 1], [1, 1, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0]]

    # test that possible_inputs_loss([1,0,1], 3) throws an ArgumentError

    @test_throws ArgumentError possible_inputs_loss([1,0,1], 3)


    @test input_probability([1,1,0,1], [1,1,0,1], QuantumDot(efficiency=1)) == 1
    @test input_probability([1,1,0,1], [1,0,0,1], QuantumDot(efficiency=1)) == 0
    @test input_probability([1,1,0,1], [1,1,0,1], QuantumDot(efficiency=0)) == 0

end