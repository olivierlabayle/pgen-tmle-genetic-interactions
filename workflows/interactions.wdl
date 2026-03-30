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
        String batch_size = "20"

        String npcs = "10"
        String approx_pca = "false"

        String ip_values = "1000 50 0.05"
        String maf = "0.01"

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
        call ld_prune {
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
        
        File plink_bed_files = ld_prune.ld_pruned_fileset.bed
        File plink_bim_files = ld_prune.ld_pruned_fileset.bim
        File plink_fam_files = ld_prune.ld_pruned_fileset.fam

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
        # Those lines aimed to send only the relevant chromosomes to the estimate_interaction call
        # However the contains function does not seem to be implemented in cromwell so I aborted
        # call extract_chroms {
        #     input: 
        #         batch_file = interaction_batch
        # }

        # scatter (pgen_fileset in pgen_filesets) {
        #     PGENFileset keep = if (contains(extract_chroms.chroms, pgen_fileset.chr)) then pgen_fileset else None
        # }

        # Array[PGENFileset] selected_pgen_filesets = select_all(keep)

        # scatter (selected_pgen_fileset in selected_pgen_filesets) {
        #     File selected_pgen_files = selected_pgen_fileset.pgen
        #     File selected_pvar_files = selected_pgen_fileset.pvar
        #     File selected_psam_files = selected_pgen_fileset.psam
        # }
        # Estimate interactions for the batch
        
        call estimate_interaction {
            input:
                docker_image = docker_image,
                julia_cmd = get_julia_cmd.julia_cmd,
                covariates_file = covariates_file,
                pcs_file = merge_and_pca.eigenvec,
                chromosomes = chromosomes,
                pgen_files = pgen_files,
                pvar_files = pvar_files,
                psam_files = psam_files,
                interaction_batch_file = interaction_batch,
                phenotype = phenotype,
                covariates = covariates,
                confounders = confounders
        }
    }

    output {
        Array[File] interaction_batches = generate_interaction_batches.interaction_batches
        Array[PLINKFileset] ld_pruned_filesets = ld_prune.ld_pruned_fileset
        File eigenvec = merge_and_pca.eigenvec
        File eigenval = merge_and_pca.eigenval
        PLINKFileset merged_fileset = merge_and_pca.merged_ld_pruned_fileset
    }
}

task estimate_interaction {
    input {
        String docker_image
        String julia_cmd
        File covariates_file
        File pcs_file
        Array[String] chromosomes
        Array[File] pgen_files
        Array[File] pvar_files
        Array[File] psam_files
        File interaction_batch_file
        String phenotype
        Array[String] covariates
        Array[String] confounders
    }

    command <<<
        for f in ~{sep=" " pgen_files}; do
            echo "${f%.bed}"
        done > pgen_files.txt

        for f in ~{sep=" " chromosomes}; do
            echo "${f%.bed}"
        done > chromosomes.txt

        covariate_opt="~{sep="," covariates}"
        if [ -z "${covariate_opt}" ]; then
            covariate_opt=" --covariates ${covariate_opt}"
        fi

        confounders_opt="~{sep="," confounders}"
        if [ -z "${confounders_opt}" ]; then
            confounders_opt=" --confounders $confounders_opt}"
        fi

        ~{julia_cmd} estimate-interactions \
            pgen_files.txt \
            chromosomes.txt \
            ~{covariates_file} \
            ~{pcs_file} \
            ~{interaction_batch_file} \
            --phenotype ~{phenotype}${covariate_opt}${confounders_opt}
    >>>

    output {
        File interaction_estimates = ""
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


task ld_prune {
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
    >>>

    output {
        PLINKFileset ld_pruned_fileset = object {
            chr: "${chr}",
            bed: "ld_pruned.no_proximal.chr_${chr}.bed",
            bim: "ld_pruned.no_proximal.chr_${chr}.bim",
            fam: "ld_pruned.no_proximal.chr_${chr}.fam"
        }
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