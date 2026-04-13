module TestEstimate

using Test
using PgenInteractions
using JSON

TESTDIR = joinpath(pkgdir(PgenInteractions), "test")

@testset "Test estimate_interactions" begin
    tmpdir = mktempdir()

    variants_list = joinpath(tmpdir, "variants_list.txt")
    open(variants_list, "w") do io
        println(io, joinpath(TESTDIR, "assets", "raw", "variants_of_interest.chr_1.raw"))
        println(io, joinpath(TESTDIR, "assets", "raw", "variants_of_interest.chr_3.raw"))
    end

    phenotype = "SEVERE_PNEUMONIA"
    output_prefix = joinpath(tmpdir, "estimates")

    copy!(ARGS, [
        "estimate-interactions",
        variants_list,
        joinpath(TESTDIR, "assets", "covariates", "covariates.csv"),
        joinpath(TESTDIR, "assets", "pcs", "ld_pruned.no_proximal.all_chr.eigenvec"),
        joinpath(TESTDIR, "assets", "interactions.batch_1.tsv"),
        "--phenotype", phenotype,
        "--covariates", "AGE",
        "--confounders", "SEX",
        "--positivity-constraint", "0",
        "--estimator-config", "ose",
        "--output-prefix", output_prefix
    ])

    julia_main()

    json_results = JSON.parsefile(string(output_prefix, ".json"))
    json_result_1 = first(json_results)
    treatment_linked_variant = Dict(
        "chr1:183905563:G:A" => ["chr1:40310265:G:A"],
        "chr3:171645187:A:G" => [],
        "chr1:40310265:G:A" => ["chr1:183905563:G:A"]
    )
    for json_result in json_results
        for estimand in json_result["OSE_GLMNET_GLMNET"]["estimand"]["args"]
            @test estimand["outcome"] == phenotype
            @test estimand["outcome_extra_covariates"] == ["AGE"]
            for (treatment, confounders) in estimand["treatment_confounders"]
                @test confounders == sort(vcat(["PC$i" for i in 1:10], treatment_linked_variant[treatment], "SEX"))
            end
        end
    end
end

end

true