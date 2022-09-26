@testset "partition types" begin

    n = 2
    m = 2
    @test is_compatible(ModeList([1,2], m),ModeList([1,2],m))
    @test is_compatible(ModeList([1,2], m),ModeList([2,1],m))
    @test_throws ArgumentError is_compatible(ModeList([1,2], m),ModeList([2,1])) == ArgumentError


end
