# pgen-tmle-genetic-interactions

Run Semi-parametric estimation of pairwise interactions from PGEN files.

## Running the Workflow on the UKB Rap

### Workflow Inputs

- `docker_image` (default: "olivierlabayle/pgen-tmle-interactions:main"): The docker image used to execute most tasks.
- `julia_use_sysimage` (default: "true"): For faster Julia execution (only disable for development.) 
- `julia_threads` (default: "auto"): Number of threads to use in Julia processes.
- `phenotype`: The phenotype under investigation, mut have a matching column in the `covariates_file`.
- `covariates` (default: ["AGE", "SEX"]): Covariates to be used by the predictive models in the TMLE step.
- `confounders` (default: []): Aditional confounders to add to the estimation process beside principal components.
- `covariates_file`: A tab separated file containing covariates and phenotype.
- `pgen_filesets`: One per chromosome, see `config/example_config.json` for an example.
- `variants_file`: A 3 columns tab separated file of variants for which all pairwise interactionswill be estimated. The columns must be:
  - `CHROM`: The chromosome without a "chr" prefix, e.g., 1, 2, 3
  - `POS`: The position of the variant.
  - `ID`: The variant unique ID.
- `batch_size` (default: "20"): For faster execution, estimation is parallelized over batches of the given size. 
- `npcs` (default: "10"): Number of principal components to use for population stratification adjustment. 
- `approx_pca` (default "true"): See [plink2's documentation](https://www.cog-genomics.org/plink/2.0/strat).
- `ip_values` (default: "1000 50 0.05"): Parameter used for LD pruning (seed [plink2's documentation](https://www.cog-genomics.org/plink/2.0/ld)).
- `maf` (default: "0.01"): Minor allele frequency used to filter genotypes for LD pruning.
- `positivity_constraint` (default: "0.01"): Similar to minor allele frequency but used at the genotype level to discard genetic effects with too low genotype frequency.
- `estimator_config` (default: "wtmle"): The semi-parametric estimator's configuration. See [this page](https://targene.github.io/TMLECLI.jl/stable/tmle_estimation/#Specifying-Estimators)

### Compile the Workflow

Using the `config/example_config.json` to be tailored to your need:

```
java -jar $DX_COMPILER_PATH compile workflows/interactions.wdl \
    -f -project $RAP_PROJECT_ID \
    -extras config/extras.json \
    -reorg \
    -folder /workflows/interactions \
    -inputs config/example_config.json
```

where:

- `DX_COMPILER_PATH` is set as per the installation's instruction.
- `RAP_PROJECT_ID` is your UK Biobank RAP project ID.

### Run the Workflow

Outputs will be generated in `/interactions_outputs`:

```
dx run -y \
-f config/example_config.dx.json \
--priority high \
--destination /interactions_outputs/ \
--preserve-job-outputs \
/workflows/interactions/interactions
```

## Building the Docker Image

```bash
docker build --platform linux/amd64 -t olivierlabayle/pgen-tmle-interactions:main -f docker/Dockerfile .
```

## Running the Docker Image


```bash
docker run -it --platform linux/amd64 olivierlabayle/pgen-tmle-interactions:main /bin/bash
```