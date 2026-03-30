module TestE2E

config = if Sys.isapple()
    "-Dconfig.file=config/cromwell.macOS-dev.conf"
else
    "-Dconfig.file=config/cromwell.local.conf"
end

cmd = Cmd([
    "java", config,
    "-jar", ENV["CROMWELL_PATH"],
    "run", joinpath(PKGDIR, "workflows", "interactions.wdl"),
    "--inputs", joinpath(TESTDIR, "assets", "config", "gwas.bygroup.json"),
    "--options", joinpath(TESTDIR, "assets", "config", "gwas.bygroup.options.json")

])

# Run the workflow from the package directory
cd(PKGDIR) do
    run(cmd)
end

end

true