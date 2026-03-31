version 1.0

import "structs.wdl"

workflow interactions {
    input {
        String docker_image = "olivierlabayle/pgen-tmle-interactions:main"
        String julia_use_sysimage = "true"
        String julia_threads = "auto"

        String phenotype
        Array[String] covariates = ["AGE", "SEX"]
        Array[String] confounders = []
        File covariates_file
        Array[PGENFileset] pgen_filesets

        File variants_file
        String batch_size = "200"

        String npcs = "10"
        String approx_pca = "true"

        String ip_values = "1000 50 0.05"
        String maf = "0.01"

        String positivity_constraint = "0.01"
        String estimator_config = "wtmle"

    }

    call get_julia_cmd as get_julia_cmd {
        input:
            use_sysimage = julia_use_sysimage,
            threads = julia_threads
    }

    call generate_interaction_batches {
        input:
            docker_image = docker_image,
            julia_cmd = get_julia_cmd.julia_cmd,
            variants_file = variants_file,
            batch_size = batch_size
    }

    scatter (pgen_fileset in pgen_filesets) {
        call ld_prune_and_extract_variants_of_interest {
            input:
                docker_image = docker_image,
                chr = pgen_fileset.chr,
                pgen_file = pgen_fileset.pgen,
                pvar_file = pgen_fileset.pvar,
                psam_file = pgen_fileset.psam,
                variants_file = variants_file,
                ip_values = ip_values,
                maf = maf
        }
        
        File plink_bed_files = ld_prune_and_extract_variants_of_interest.ld_pruned_fileset.bed
        File plink_bim_files = ld_prune_and_extract_variants_of_interest.ld_pruned_fileset.bim
        File plink_fam_files = ld_prune_and_extract_variants_of_interest.ld_pruned_fileset.fam

        File pgen_files = pgen_fileset.pgen
        File pvar_files = pgen_fileset.pvar
        File psam_files = pgen_fileset.psam
        String chromosomes = pgen_fileset.chr
    }

    call merge_and_pca {
        input:
            docker_image = docker_image,
            plink_bed_files = plink_bed_files,
            plink_bim_files = plink_bim_files,
            plink_fam_files = plink_fam_files,
            npcs = npcs,
            approx_pca = approx_pca
    }

    scatter (interaction_batch in generate_interaction_batches.interaction_batches) {
        call estimate_interaction {
            input:
                docker_image = docker_image,
                julia_cmd = get_julia_cmd.julia_cmd,
                covariates_file = covariates_file,
                pcs_file = merge_and_pca.eigenvec,
                variant_data = select_all(ld_prune_and_extract_variants_of_interest.variants_of_interest),
                interaction_batch_file = interaction_batch,
                phenotype = phenotype,
                covariates = covariates,
                confounders = confounders,
                positivity_constraint = positivity_constraint,
                estimator_config = estimator_config
        }
    }

    call generate_outputs {
        input:
            docker_image = docker_image,
            julia_cmd = get_julia_cmd.julia_cmd,
            hdf5_files = select_all(estimate_interaction.hdf5_interaction_estimates)
    }

    output {
        Array[File] interaction_batches = generate_interaction_batches.interaction_batches
        Array[PLINKFileset] ld_pruned_filesets = ld_prune_and_extract_variants_of_interest.ld_pruned_fileset
        File eigenvec = merge_and_pca.eigenvec
        File eigenval = merge_and_pca.eigenval
        PLINKFileset merged_fileset = merge_and_pca.merged_ld_pruned_fileset
        Array[File] json_estimates = select_all(estimate_interaction.json_interaction_estimates)
        File hdf5_output = generate_outputs.hdf5_output
        File yaml_output = generate_outputs.yaml_output
        File qq_output = generate_outputs.qq_output
    }
}


task generate_outputs {
    input {
        String docker_image
        String julia_cmd
        Array[File] hdf5_files
    }

    command <<<
        for f in ~{sep=" " hdf5_files}; do
            ln -s "$f" .
        done

        ~{julia_cmd} make-outputs interactions.batch_
    >>>

    output {
        File hdf5_output = "results.hdf5"
        File yaml_output = "results.summary.yaml"
        File qq_output = "QQ.png"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

}

task estimate_interaction {
    input {
        String docker_image
        String julia_cmd
        File covariates_file
        File pcs_file
        Array[File] variant_data
        File interaction_batch_file
        String phenotype
        Array[String] covariates
        Array[String] confounders
        String positivity_constraint = "0.01"
        String estimator_config = "wtmle"
    }

    String output_prefix = basename(interaction_batch_file, ".tsv")

    command <<<
        for f in ~{sep=" " variant_data}; do
            echo "${f}"
        done > variant_files.txt

        covariate_opt="~{sep="," covariates}"
        if [ -n "${covariate_opt}" ]; then
            covariate_opt=" --covariates ${covariate_opt}"
        fi

        confounders_opt="~{sep="," confounders}"
        if [ -n "${confounders_opt}" ]; then
            confounders_opt=" --confounders $confounders_opt}"
        fi

        ~{julia_cmd} estimate-interactions \
            variant_files.txt \
            ~{covariates_file} \
            ~{pcs_file} \
            ~{interaction_batch_file} \
            --positivity-constraint ~{positivity_constraint} \
            --estimator-config ~{estimator_config} \
            --output-prefix ~{output_prefix} \
            --phenotype ~{phenotype}${covariate_opt}${confounders_opt}
    >>>

    output {
        File? json_interaction_estimates = "${output_prefix}.json"
        File? hdf5_interaction_estimates = "${output_prefix}.hdf5"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

}

task get_julia_cmd {
    input {
        String use_sysimage = "true"
        String threads = "auto"
    }
    command <<<
        julia_cmd_string="julia --project=/opt/PgenInteractions --startup-file=no"
        if [[ "~{use_sysimage}" == "true" ]]; then
            julia_cmd_string+=" --sysimage=/opt/PgenInteractions/sysimage.so"
        fi
        if [[ "~{threads}" == "auto" ]]; then
            julia_cmd_string+=" --threads=auto"
        fi
        julia_cmd_string+=" /opt/PgenInteractions/bin/run.jl"
        echo "$julia_cmd_string"
    >>>

    output {
        String julia_cmd = read_string(stdout())
    }
}

task extract_chroms {
  input {
    File batch_file
  }

  command <<<
    # extract both columns, remove 'chr', flatten, sort & unique
    tail -n +2 ~{batch_file} | cut -f1 > col1.txt
    tail -n +2 ~{batch_file} | cut -f4 > col2.txt
    cat col1.txt col2.txt | sed 's/^chr//' | sort | uniq > out.txt
  >>>

  output {
    Array[String] chroms = read_lines("out.txt")
  }
}

task merge_and_pca {
    input {
        String docker_image
        Array[File] plink_bed_files
        Array[File] plink_bim_files
        Array[File] plink_fam_files
        String npcs = "10"
        String approx_pca = "false"
    }

    command <<<
        # Merge filesets
        for f in ~{sep=" " plink_bed_files}; do
            echo "${f%.bed}"
        done > merge_list.txt

        plink \
            --biallelic-only \
            --merge-list merge_list.txt \
            --make-bed \
            --out ld_pruned.no_proximal.all_chr

        approx_option=""
        if [[ "~{approx_pca}" == "true" ]]; then
            approx_option=" approx"
        fi

        plink2 \
            --bfile ld_pruned.no_proximal.all_chr \
            --pca ~{npcs}${approx_option} \
            --out ld_pruned.no_proximal.all_chr
    >>>

    output {
        PLINKFileset merged_ld_pruned_fileset = object {
            chr: "all",
            bed: "ld_pruned.no_proximal.all_chr.bed",
            bim: "ld_pruned.no_proximal.all_chr.bim",
            fam: "ld_pruned.no_proximal.all_chr.fam"
        }

        File eigenvec = "ld_pruned.no_proximal.all_chr.eigenvec"
        File eigenval = "ld_pruned.no_proximal.all_chr.eigenval"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}


task ld_prune_and_extract_variants_of_interest {
    input {
        String docker_image
        String chr
        File pgen_file
        File pvar_file
        File psam_file
        File variants_file
        String ip_values = "1000 50 0.05"
        String maf = "0.01"
    }

    command <<<
        pgen_prefix=$(dirname "~{pgen_file}")/$(basename "~{pgen_file}" .pgen)

        plink2 \
            --pfile ${pgen_prefix} \
            --rm-dup force-first \
            --indep-pairwise ~{ip_values}
        # Always exclude high LD regions stored in the docker image at /opt/PgenInteractions/exclude_b38.txt
        plink2 \
            --pfile ${pgen_prefix} \
            --extract plink2.prune.in \
            --maf ~{maf} \
            --make-bed \
            --exclude range /opt/PgenInteractions/assets/exclude_b38.txt \
            --out ld_pruned.chr_~{chr}
        # Exlude variants around queried variants to avoid proximal contamination
        awk 'BEGIN{OFS="\t"}
            NR==1 {
            for (i=1; i<=NF; i++) {
                if ($i=="CHROM") c=i
                else if ($i=="POS") p=i
                else if ($i=="ID") id=i
            }
            next
            }
            {
            start = $p - 500000
            if (start < 0) start = 0
            end = $p + 500000
            print $c, start, end
        }' ~{variants_file} > exclude_ranges.txt

        plink2 \
            --bfile ld_pruned.chr_~{chr} \
            --exclude range exclude_ranges.txt \
            --make-bed \
            --out ld_pruned.no_proximal.chr_~{chr}

        # Extract raw data for the variants of interest
        awk -F'\t' '
            NR==1 {
            for (i=1; i<=NF; i++) if ($i=="ID") col=i
            next
            }
            { print $col }
        ' ~{variants_file} > variants_of_interest.txt

        plink2 \
            --pfile ${pgen_prefix} \
            --extract variants_of_interest.txt \
            --export A include-alt \
            --out variants_of_interest.chr_~{chr} || echo "No variant of interest on this chromosome."
    >>>

    output {
        PLINKFileset ld_pruned_fileset = object {
            chr: "${chr}",
            bed: "ld_pruned.no_proximal.chr_${chr}.bed",
            bim: "ld_pruned.no_proximal.chr_${chr}.bim",
            fam: "ld_pruned.no_proximal.chr_${chr}.fam"
        }
        File? variants_of_interest = "variants_of_interest.chr_${chr}.raw"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}

task generate_interaction_batches {
    input {
        String docker_image
        String julia_cmd
        File variants_file
        String batch_size = "20"
    }

    command <<<
        ~{julia_cmd} generate-interactions \
            ~{variants_file} \
            --output-prefix interactions \
            --batch-size ~{batch_size}
    >>>

    output {
        Array[File] interaction_batches = glob("interactions.*.tsv")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}