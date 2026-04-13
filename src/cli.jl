
function cli_settings()
    s = ArgParseSettings(
        description="Pgen-Interactions",
        add_version = true,
        commands_are_required = false,
        version=string(pkgversion(PgenInteractions))
    )

    @add_arg_table! s begin
        "generate-interactions"
            action = :command
            help = "Generate interaction files to be estimated in parallel."

        "estimate-interactions"
            action = :command
            help = "Estimate all interactions in the batch file"

        "make-outputs"
            action = :command
            help = "Generate pipeline outputs."

    end

    @add_arg_table! s["generate-interactions"] begin
        "variants-file"
            arg_type = String
            required = true
            help = "Path to variants file containing 3 columns: CHROM, POS, ID"
    
        "--output-prefix"
            arg_type = String
            required = false
            help = "Output prefix of interaction batches."
            default = "interactions"

        "--batch-size"
            arg_type = Int
            required = false
            help = "Maximum number of interactions to be estimaed per batch."
            default = 20

        "--adjustment-window-kb"
            arg_type = Float64
            required = false
            help = "All variants in the window will be included in the confouders set."
            default = 1000
    end

    @add_arg_table! s["estimate-interactions"] begin
        "variant-files"
            arg_type = String
            help = "Path to list of variant files"

        "covariates-file"
            arg_type = String
            help = "Path to covariates"

        "pcs-file"
            arg_type = String
            help = "Path to pcs"

        "interaction-batch-file"
            arg_type = String
            help = "Path to interaction batch to estimate"
        
        "--positivity-constraint"
            arg_type = Float64
            default = 0.01
            help = "Positivity constraint to apply on estimands"

        "--phenotype"
            arg_type = String
            default = "Y"
            help = "Target phenotype"

        "--covariates"
            arg_type = String
            default = nothing
            help = "Comma separated list of extra covariates"
        
        "--confounders"
            arg_type = String
            default = nothing
            help = "Comma separated list of extra confounders"

        "--estimator-config"
            arg_type = String
            default = nothing
            help = "Estimator's configuration"

        "--output-prefix"
            arg_type = String
            default = "estimates"
            help = "Output prefix"
    end

    @add_arg_table! s["make-outputs"] begin
        "input-prefix"
            help = "Input prefix to the results file."
            required = true

        "--output-prefix"
            help = "Prefix to output plots."
            default = "."

        "--verbosity"
            help = "Logging level."
            arg_type = Int
            default = 0
    end

    return s
end

function julia_main()::Cint
    settings = parse_args(ARGS, cli_settings())
    cmd = settings["%COMMAND%"]
    @info "Running Pgen-Interactions CLI: $cmd"
    cmd_settings = settings[cmd]
    if cmd == "generate-interactions"
        generate_interaction_batches(
            cmd_settings["variants-file"];
            output_prefix=cmd_settings["output-prefix"],
            batch_size=cmd_settings["batch-size"],
            adjustment_window_kb=cmd_settings["adjustment-window-kb"]
        )
    elseif cmd == "estimate-interactions"
        estimate_interactions(
            cmd_settings["variant-files"],
            cmd_settings["covariates-file"],
            cmd_settings["pcs-file"],
            cmd_settings["interaction-batch-file"];
            phenotype=cmd_settings["phenotype"],
            covariates=cmd_settings["covariates"],
            confounders=cmd_settings["confounders"],
            positivity_constraint=cmd_settings["positivity-constraint"],
            estimator_config=cmd_settings["estimator-config"],
            output_prefix=cmd_settings["output-prefix"]
        )
    elseif cmd == "make-outputs"
        make_outputs(
            cmd_settings["input-prefix"];
            output_prefix=cmd_settings["output-prefix"],
            verbosity=cmd_settings["verbosity"],
        )
    else
        throw(ArgumentError(string("Unknown command: ", cmd)))
    end
    return 0
end