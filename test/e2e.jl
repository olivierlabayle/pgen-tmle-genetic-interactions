module TestE2E

using PgenInteractions
using Test

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

end

true