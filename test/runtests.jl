using Test

@testset "All Tests" begin
    @test include("estimands.jl")
    @test include("estimate.jl")
    @test include("e2e.jl")
end