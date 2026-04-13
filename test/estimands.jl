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
        "--batch-size", string(batchsize),
        "--adjustment-window-kb", string("500")
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
    @test nrow(unique(all_interactions, [:ID_1, :ID_2])) == nrow(all_interactions)
    @test names(all_interactions) == [
        "CHROM_1",
        "POS_1",
        "ID_1",
        "LINKED_1",
        "CHROM_2",
        "POS_2",
        "ID_2",
        "LINKED_2"
    ]
    example_variant = "chr1:155072528:T:C"
    linked_data = subset(all_interactions, :ID_1 => x -> x .== example_variant).LINKED_1
    linked_variants = only(unique(split.(linked_data, ";&;")))
    @test example_variant ∉ linked_variants
    @test Set([
        "chr1:155093927:C:T",
        "chr1:155177242:G:A",
        "chr1:155192003:C:T",
        "chr1:155197995:A:G",
        "chr1:155212198:G:A"]) == Set(linked_variants)
end

end

true