module TestEstimate

using Test
using PgenInteractions
using JSON

TESTDIR = joinpath(pkgdir(PgenInteractions), "test")

@testset "Test estimate_interactions" begin
    tmpdir = mktempdir()

    pgen_list = joinpath(tmpdir, "pgen.txt")
    open(pgen_list, "w") do io
        println(io, joinpath(TESTDIR, "assets", "pgen", "chr1.qced.pgen"))
        println(io, joinpath(TESTDIR, "assets", "pgen", "chr3.qced.pgen"))
    end
    chrom_list = joinpath(tmpdir, "chrom.txt")
    open(chrom_list, "w") do io
        println(io, "1")
        println(io, "3")
    end

    phenotype = "SEVERE_PNEUMONIA"
    output_prefix = joinpath(tmpdir, "estimates")

    copy!(ARGS, [
        "estimate-interactions",
        pgen_list,
        chrom_list,
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
    for json_result in json_results
        for estimand in json_result["OSE_GLMNET_GLMNET"]["estimand"]["args"]
            @test estimand["outcome"] == phenotype
            @test estimand["outcome_extra_covariates"] == ["AGE"]
            for (treatment, confounders) in estimand["treatment_confounders"]
                @test confounders == sort(vcat(["PC$i" for i in 1:10], "SEX"))
            end
        end
    end
end

end

true