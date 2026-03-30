version 1.0

workflow interactions {
    input {
        String docker_image = "olivierlabayle/pgen-tmle-interactions:main"
        File variants_file
        Int batch_size = 20
        String use_sysimage = "true"
        String threads = "auto"
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
        julia_cmd_string+=" /opt/PgenInteractions/run.jl"
        echo "$julia_cmd_string"
    >>>

    output {
        String julia_cmd = read_string(stdout())
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