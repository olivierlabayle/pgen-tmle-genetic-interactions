using Test

@testset "All Tests" begin
    @test include("estimands.jl")
end