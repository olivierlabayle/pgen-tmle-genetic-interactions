using Test
using PgenInteractions

TESTDIR = joinpath(pkgdir(PgenInteractions), "test")

@testset "All Tests" begin
    @test include(joinpath(TESTDIR, "estimands.jl"))
    @test include(joinpath(TESTDIR, "estimate.jl"))
    @test include(joinpath(TESTDIR, "e2e.jl"))
end