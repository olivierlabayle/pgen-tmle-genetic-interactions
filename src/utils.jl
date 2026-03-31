function files_matching_prefix(prefix)
    directory, _prefix = splitdir(prefix)
    _directory = directory == "" ? "." : directory

    return map(
        f -> joinpath(directory, f),
        filter(
            f -> startswith(f, _prefix), 
            readdir(_directory)
        )
    )
end

function make_filepath_from_prefix(prefix; filename="dataset.arrow")
    return if isdir(prefix)
        joinpath(prefix, filename)
    else
        dir, _prefix = splitdir(prefix)
        joinpath(dir, string(_prefix, filename))
    end
end

###############################################################################
###                    P-VALUES WITH ERROR MANAGEMENT                       ###
###############################################################################

function pvalue_or_nan(Ψ̂)
    return try
        pvalue(significance_test(Ψ̂))
    catch
        NaN
    end
end

function pvalue_or_nan(Ψ̂, Ψ₀)
    return try
        pvalue(significance_test(Ψ̂, Ψ₀))
    catch
        NaN
    end
end

pvalue_or_nan(Ψ̂::TMLECLI.FailedEstimate, Ψ₀=nothing) = NaN