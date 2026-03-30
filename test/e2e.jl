module TestE2E

using PgenInteractions
using Test
using CSV
using DataFrames

PKGDIR = pkgdir(PgenInteractions)

TESTDIR = joinpath(PKGDIR, "test")

config = if Sys.isapple()
    "-Dconfig.file=config/cromwell.macOS-dev.conf"
else
    "-Dconfig.file=config/cromwell.local.conf"
end

cmd = Cmd([
    "java", config,
    "-jar", ENV["CROMWELL_PATH"],
    "run", joinpath(PKGDIR, "workflows", "interactions.wdl"),
    "--inputs", joinpath(TESTDIR, "config", "config.json"),
    "--options", joinpath(TESTDIR, "config", "config.options.json")
])

# Run the workflow from the package directory
cd(PKGDIR) do
    run(cmd)
end

results_dirs = readdir(joinpath(PKGDIR, "interactions_workflow_outputs", "interactions"), join=true)
results_dir = results_dirs[argmax(mtime(d) for d in results_dirs)]

# Check batches generation
batch_files_execution_dir = joinpath(results_dir, "call-generate_interaction_batches/execution")
batch_files_output_dir = only(readdir(batch_files_execution_dir, join=true))
batch_files = readdir(batch_files_output_dir, join=true)
interactions_df = mapreduce(batch_file -> CSV.read(batch_file, DataFrame), vcat, batch_files)
@test nrow(interactions_df) == 3

# Test LD pruning removes variants in interactions
ld_pruning_dir = joinpath(results_dir, "call-ld_prune")
for shard in 0:2
    ld_output_dir = joinpath(ld_pruning_dir, "shard-$shard", "execution")
    bim = CSV.read(
        joinpath(ld_output_dir, "ld_pruned.no_proximal.chr_$(shard+1).bim"), 
        DataFrame; 
        header=[:CHROM, :ID, :CM, :POS, :REF, :ALT]
    )
    @test length(intersect(bim.ID, interactions_df.ID_1)) == 0
    @test length(intersect(bim.ID, interactions_df.ID_2)) == 0
end

# Test merging and PCA
pca_output_dir = joinpath(results_dir, "call-merge_and_pca", "execution")
pruned_variants_df = CSV.read(
    joinpath(pca_output_dir, "ld_pruned.no_proximal.all_chr.bim"),
    DataFrame;
    header=[:CHROM, :ID, :CM, :POS, :REF, :ALT]
)
@test nrow(pruned_variants_df) > 15
@test length(readlines(joinpath(pca_output_dir, "ld_pruned.no_proximal.all_chr.eigenval"))) == 10
pcs = CSV.read(joinpath(pca_output_dir, "ld_pruned.no_proximal.all_chr.eigenvec"), DataFrame)
@test names(pcs) == ["#FID", "IID", ("PC$i" for i in 1:10)...]



end

true