function generate_interaction_batches(variants_file; output_prefix="interactions", batch_size=20)
    # Read variants
    variants = CSV.read(variants_file, DataFrame;select=[:CHROM, :POS, :ID])
    # Generate Interactions
    variants = NamedTuple.(eachrow(variants))
    n_variants = length(variants)
    variants_interactions = NamedTuple[]
    for variant_idx_1 in 1:n_variants
        for variant_idx_2 in (variant_idx_1+1):n_variants
            variant_1 = variants[variant_idx_1]
            variant_2 = variants[variant_idx_2]
            push!(variants_interactions, (
                CHROM_1=variant_1.CHROM, 
                POS_1=variant_1.POS, 
                ID_1=variant_1.ID,
                CHROM_2=variant_2.CHROM, 
                POS_2=variant_2.POS, 
                ID_2=variant_2.ID,
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