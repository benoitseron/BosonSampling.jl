@testset "sampling structures" begin
    
    params = PartitionSamplingParameters(n = 10, m = 10)

    # set_input!(params)
    set_interferometer!(build_loop(LoopSamplingParameters(m=10)), params)

    compute_probability!(params)
end