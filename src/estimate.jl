function update_chr_to_variant!(chr_to_variants, key, val)
    if haskey(chr_to_variants, key)
        push!(chr_to_variants[key], val)
    else
        chr_to_variants[key] = Set([val])
    end
end

function get_chr_to_variants(interaction_variants)
    chr_to_variants = Dict()
    for row in eachrow(interaction_variants)
        update_chr_to_variant!(chr_to_variants, row.CHROM_1, row.ID_1)
        update_chr_to_variant!(chr_to_variants, row.CHROM_2, row.ID_2)
    end
    return chr_to_variants
end

function extract_variants_from_pgens(interaction_variants, chr_to_pgen)
    chr_to_variants = get_chr_to_variants(interaction_variants)
    
    # Extract all relevant variants to a table using plink2
    variants_dfs = []
    for (chr, variants) in chr_to_variants
        pgen_prefix = replace(chr_to_pgen[string(chr)], ".pgen" => "")
        tmpdir = mktempdir()
        extract_list = joinpath(tmpdir, "variants.txt")
        open(extract_list, "w") do io
            for v in variants
                println(io, v)
            end
        end
        outprefix = joinpath(tmpdir, string("chr", chr))
        run(`plink2 \
            --pfile $pgen_prefix \
            --extract $extract_list \
            --export A include-alt \
            --out $outprefix
        `)
        push!(
            variants_dfs,
            CSV.read(string(outprefix, ".raw"), DataFrame; missingstring="NA", drop=[:PAT, :MAT, :SEX, :PHENOTYPE])
        )
    end
    variants_df = popfirst!(variants_dfs)
    while length(variants_dfs) > 0
        variants_df = innerjoin(
            variants_df, 
            popfirst!(variants_dfs),
            on=[:FID, :IID]
        )
    end
    # Rename variant Ids to their original names
    for name in names(variants_df)
        if name ∉ ["FID", "IID"]
            original_name = join(split(name, "_")[1:end-1], "_")
            rename!(variants_df, name => original_name)
        end
    end
    return variants_df
end

make_list_from_option(opt::Nothing) = []

make_list_from_option(x) = split(x, ",")

function make_dataset(interaction_variants, chr_to_pgen, pcs_file, covariates_file; 
    covariates=[], 
    confounders=[], 
    phenotype="Y"
    )
    dataset = extract_variants_from_pgens(interaction_variants, chr_to_pgen)
    dataset = innerjoin(
        dataset, 
        CSV.read(pcs_file, DataFrame), 
        on = ["FID" => "#FID", "IID"]
    )
    return innerjoin(
        dataset, 
        CSV.read(
            covariates_file, 
            DataFrame; 
            select=["FID", "IID", covariates..., confounders..., phenotype]
        ),
        on=[:FID, :IID]
    )
end

function estimate_interactions(
    pgen_file_list,
    chrom_list,
    covariates_file, 
    pcs_file,
    interaction_batch_file;
    phenotype="Y",
    covariates=nothing,
    confounders=nothing
    )
    positivity_constraint = 0.01
    estimators_config = "wtmle"
    output_prefix = replace(basename(interaction_batch_file), ".tsv" => "")

    phenotype = "SEVERE_COVID_19"
    covariates = "AGE,SEX"
    confounders = nothing
    pgen_file_list = "cromwell-executions/interactions/cec5e018-21e7-419f-97ff-575e5cd89e74/call-estimate_interaction/shard-1/execution/pgen_files.txt"
    chrom_list = "cromwell-executions/interactions/cec5e018-21e7-419f-97ff-575e5cd89e74/call-estimate_interaction/shard-1/execution/chromosomes.txt"
    interaction_batch_file = "cromwell-executions/interactions/cec5e018-21e7-419f-97ff-575e5cd89e74/call-estimate_interaction/shard-1/inputs/-891516308/interactions.batch_2.tsv"
    covariates_file = "cromwell-executions/interactions/cec5e018-21e7-419f-97ff-575e5cd89e74/call-estimate_interaction/shard-1/inputs/250885478/covariates.csv"
    pcs_file = "cromwell-executions/interactions/cec5e018-21e7-419f-97ff-575e5cd89e74/call-estimate_interaction/shard-1/inputs/390391270/ld_pruned.no_proximal.all_chr.eigenvec"
    
    confounders = make_list_from_option(confounders)
    covariates = make_list_from_option(covariates)
    chr_to_pgen = Dict(chr => file for (chr, file) in zip(readlines(chrom_list), readlines(pgen_file_list)))
    interaction_variants = CSV.read(interaction_batch_file, DataFrame)
    dataset = make_dataset(interaction_variants, chr_to_pgen, pcs_file, covariates_file; 
        covariates=covariates, 
        confounders=confounders, 
        phenotype=phenotype
    )
    pcs = filter(startswith("PC"), names(dataset))
    all_confounders = vcat(confounders, pcs)
    estimands = mapreduce(vcat, eachrow(interaction_variants)) do row
        factorialEstimand(
            AIE,
            (row.ID_1, row.ID_2),
            phenotype;
            confounders=all_confounders,
            positivity_constraint=positivity_constraint,
            outcome_extra_covariates=covariates,
            dataset=dataset
        )
    end

    runner = Runner(dataset;
        estimands_config=TMLE.Configuration(estimands=estimands), 
        estimators_spec=estimators_config, 
        outputs=TMLECLI.Outputs(json=string(output_prefix, ".json"), hdf5=string(output_prefix, ".json")), 
    )
    runner()

    return 0
end