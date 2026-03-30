module TestEstimands

using Test
using PgenInteractions
using CSV
using DataFrames

TESTDIR = joinpath(pkgdir(PgenInteractions), "test")

@testset "Test generate_estimands" begin
    variants_file = joinpath(TESTDIR, "assets", "clumps.tsv")
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "interactions")
    outcome = "SEVERE_COVID_19"
    batchsize = 20
    copy!(ARGS,
        [
        "generate-interactions",
        variants_file,
        "--output-prefix", output_prefix,
        "--batch-size", string(batchsize)
        ]
    )
    julia_main()
    batch_files = readdir(tmpdir; join=true)
    all_interactions = mapreduce(vcat, batch_files) do batch_file
        interactions = CSV.read(batch_file, DataFrame)
        @test nrow(interactions) <= 20

        interactions
    end
    n_variants = countlines(variants_file) - 1
    @test nrow(all_interactions) == 12 * 11 / 2
    @test unique(all_interactions, [:ID_1, :ID_2]) == all_interactions
    @test names(all_interactions) == [
        "CHROM_1",
        "POS_1",
        "ID_1",
        "CHROM_2",
        "POS_2",
        "ID_2"
    ]
end

end

true