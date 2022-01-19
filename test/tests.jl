@testset "permanent_ryser" begin

    theoretical_permanents_fourier_matrix = [1,0,-3,0,-5,0,-105,0,81,0,6765,0,175747,0,30375,0,25219857,0,142901109,0,4548104883,0]
    #from "on the permanent of schur matrix", graham
    for n = 1:10
        U = ones(Int64,n,n)
        @test permanent_ryser(U) == factorial(n)
        U = fourier_matrix(n, normalized = false)
        @test permanent_ryser(U) ≈ theoretical_permanents_fourier_matrix[n] atol=0.00001
    end

end
