
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
            batch_size=cmd_settings["batch-size"]
        )
    else
        throw(ArgumentError(string("Unknown command: ", cmd)))
    end
    return 0
end