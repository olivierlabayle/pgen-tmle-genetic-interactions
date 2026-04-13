function get_linked_variants(variant, variants_df;adjustment_window_kb=1000)
    linked_variants = subset(
        variants_df,
        :CHROM => x -> x .== variant.CHROM,
        :POS => x -> variant.POS - adjustment_window_kb*1000 .<= x .<= variant.POS + adjustment_window_kb*1000
    ).ID
    return filter(!=(variant.ID), linked_variants)
end

function generate_interaction_batches(variants_file; output_prefix="interactions", batch_size=20, adjustment_window_kb=1000)
    # Read variants
    variants_df = CSV.read(variants_file, DataFrame; select=[:CHROM, :POS, :ID])
    # Generate Interactions
    variants = NamedTuple.(eachrow(variants_df))
    n_variants = length(variants)
    variants_interactions = NamedTuple[]
    for variant_idx_1 in 1:n_variants
        for variant_idx_2 in (variant_idx_1+1):n_variants
            variant_1 = variants[variant_idx_1]
            linked_variants_1 = get_linked_variants(variant_1, variants_df; adjustment_window_kb=adjustment_window_kb)
            variant_2 = variants[variant_idx_2]
            linked_variants_2 = get_linked_variants(variant_2, variants_df; adjustment_window_kb=adjustment_window_kb)
            push!(variants_interactions, (
                CHROM_1=variant_1.CHROM, 
                POS_1=variant_1.POS, 
                ID_1=variant_1.ID,
                LINKED_1=join(linked_variants_1, ";&;"),
                CHROM_2=variant_2.CHROM, 
                POS_2=variant_2.POS, 
                ID_2=variant_2.ID,
                LINKED_2=join(linked_variants_2, ";&;")
                )
            )
        end
    end
    # Write batches
    variants_interactions = DataFrame(variants_interactions)
    batch_idx = 1
    for chr_group in groupby(variants_interactions, [:CHROM_1, :CHROM_2])
        sort!(chr_group, [:ID_1, :ID_2])
        for batch in Iterators.partition(1:nrow(chr_group), batch_size)
            CSV.write(string(output_prefix, ".batch_", batch_idx, ".tsv"), chr_group[batch, :]; delim="\t")
            batch_idx += 1
        end
    end

    return 0
end