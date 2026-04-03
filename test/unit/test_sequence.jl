@testset "Sequence rejects non-positive dt" begin
    @test_throws ArgumentError Sequence(0.0)
    @test_throws ArgumentError Sequence(-1e-9)
end
