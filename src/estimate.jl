function hardcall!(hard_calls, dosages, index; threshold=0.2)
    dosage = dosages[index]
    ismissing(dosage) && return missing
    hard_call = round(Int, dosage)
    if abs(hard_call - dosage) < threshold
        hard_calls[index] = hard_call
    end
end

function hardcall(dosages; threshold=0.2)
    hard_calls = Vector{Union{Missing, Int}}(missing, size(dosages, 1))
    for index in eachindex(dosages)
        hardcall!(hard_calls, dosages, index; threshold=threshold)
    end
    return hard_calls
end

function extract_variants_from_raw_files(interaction_variants, variant_files; hard_call_threshold=0.2)
    linked_1 = vcat(split.(skipmissing(interaction_variants.LINKED_1), ";&;")...)
    linked_2 = vcat(split.(skipmissing(interaction_variants.LINKED_2), ";&;")...)
    variants_of_interest = unique(vcat(interaction_variants.ID_1, interaction_variants.ID_2, linked_1, linked_2))

    variants_dfs = []
    for raw_file in readlines(variant_files)
        variant_names_in_file = names(CSV.read(raw_file, DataFrame; limit=0, drop=[:FID, :IID, :PAT, :MAT, :SEX, :PHENOTYPE]))
        source_variant_names_in_file = [join(split(name, "_")[1:end-1], "_") for name in variant_names_in_file]
        variants_of_interest_in_file = intersect(variants_of_interest, source_variant_names_in_file)
        if !isempty(variants_of_interest_in_file)
            file_to_src_variant_map = Dict(in_file => src for (in_file, src) in zip(variant_names_in_file, source_variant_names_in_file) if src in variants_of_interest_in_file)
            raw_data = CSV.read(raw_file, DataFrame; missingstring="NA", select=["FID", "IID", keys(file_to_src_variant_map)...])
            rename!(raw_data, file_to_src_variant_map)
            for variant in values(file_to_src_variant_map)
                raw_data[!, variant] = hardcall(raw_data[!, variant]; threshold=hard_call_threshold)
            end
            push!(variants_dfs, raw_data)
        end
    end
    # Join all variant_data
    variants_df = popfirst!(variants_dfs)
    while length(variants_dfs) > 0
        variants_df = innerjoin(
            variants_df, 
            popfirst!(variants_dfs),
            on=[:FID, :IID]
        )
    end

    return variants_df
end

make_list_from_option(opt::Nothing) = []

make_list_from_option(x) = split(x, ",")

function make_dataset(interaction_variants, variant_files, pcs_file, covariates_file; 
    covariates=[], 
    confounders=[], 
    phenotype="Y",
    hard_call_threshold=0.2
    )
    dataset = extract_variants_from_raw_files(interaction_variants, variant_files; hard_call_threshold=hard_call_threshold)
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

confounder_string_to_list(x::Missing) = []

confounder_string_to_list(x) = split(x, ";&;")

make_variant_confounders(all_confounders, linked_variants, variants_in_interaction) = 
    Symbol.(unique(vcat(all_confounders, setdiff(confounder_string_to_list(linked_variants), variants_in_interaction))))

function generate_variants_confounders(row, all_confounders)
    variants_in_interaction = [row.ID_1, row.ID_2]
    return Dict(
        Symbol(row.ID_1) => make_variant_confounders(all_confounders, row.LINKED_1, variants_in_interaction),
        Symbol(row.ID_2) => make_variant_confounders(all_confounders, row.LINKED_2, variants_in_interaction),
    )
end

function estimate_interactions(
    variant_files,
    covariates_file, 
    pcs_file,
    interaction_batch_file;
    phenotype="Y",
    covariates=nothing,
    confounders=nothing,
    positivity_constraint=0.01,
    estimator_config="wtmle",
    output_prefix="estimates",
    hard_call_threshold=0.2
    )

    confounders = make_list_from_option(confounders)
    covariates = make_list_from_option(covariates)
    interaction_variants = CSV.read(interaction_batch_file, DataFrame)
    dataset = make_dataset(interaction_variants, variant_files, pcs_file, covariates_file; 
        covariates=covariates, 
        confounders=confounders, 
        phenotype=phenotype,
        hard_call_threshold = hard_call_threshold
    )
    pc_variables = filter(startswith("PC"), names(dataset))
    all_confounders = vcat(confounders, pc_variables)

    estimands = []
    for row in eachrow(interaction_variants)
        (string(row.ID_1) in names(dataset) && string(row.ID_2) in names(dataset)) || continue
        try 
            Ψ = factorialEstimand(
                AIE,
                (row.ID_1, row.ID_2),
                phenotype;
                confounders=generate_variants_confounders(row, all_confounders),
                positivity_constraint=positivity_constraint,
                outcome_extra_covariates=covariates,
                dataset=dataset,
                verbosity=0
            )
            push!(estimands, Ψ)
        catch err
            if err != ArgumentError("No component passed the positivity constraint.")
                throw(err)
            end
        end
    end

    if length(estimands) > 0
        runner = Runner(dataset;
            estimands_config=TMLE.Configuration(estimands=estimands), 
            estimators_spec=estimator_config, 
            outputs=TMLECLI.Outputs(json=string(output_prefix, ".json"), hdf5=string(output_prefix, ".hdf5")), 
        )
        runner()
    end

    return 0
end