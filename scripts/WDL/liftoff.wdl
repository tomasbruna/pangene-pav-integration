version 1.0

workflow liftoff {
    input {
        File targetFasta
        File referenceFasta
        File referenceGff
        Int? n_cpu = 32
        String? outputName = "lifted.gff3"
    }

    call unzip as unzipTargetFasta {
       input: in = targetFasta
    }

    call unzip as unzipReferenceFasta {
       input: in = referenceFasta
    }

    call unzip as unzipReferenceGff {
       input: in = referenceGff
    }

    call lift {
        input: targetFasta = unzipTargetFasta.unzipped,
               referenceFasta = unzipReferenceFasta.unzipped,
               referenceGff = unzipReferenceGff.unzipped,
               n_cpu=n_cpu
    }

    call processResult {
        input: lifted = lift.lifted,
               reference = unzipReferenceGff.unzipped,
               outputName = outputName
    }
    output {
        File result = processResult.result
    }
}

task unzip {
    # If a file is already unzipped, just copy it. Originally, the file
    # was only symlinked, but symlinks do not play nice with Docker
    input {
        File in
    }

    command <<<
        base=$(basename "~{in}")
        if [[ ${base##*.} == "gz" ]]; then
            gunzip "~{in}" -c > "${base%.*}"
        else
            cp "~{in}" "$base"
        fi
    >>>

    output {
       File unzipped = basename(in, ".gz")
    }
}

task lift {

    input {
        File targetFasta
        File referenceFasta
        File referenceGff
        Int? n_cpu = 64
    }

    command <<<
        set -e
        echo -e "gene\nmRNA" > types
        liftoff ~{targetFasta} ~{referenceFasta} -g ~{referenceGff} \
            -o intermediate.gff3_polished -p ~{n_cpu} -f types
    >>>


    output {
       File lifted = "intermediate.gff3_polished"
    }
    
    runtime {
        docker: "quay.io/biocontainers/liftoff@sha256:63d9a69375519259f155e2f0b0a61b4c95287684f324c7193e2cead7e4ef5894"
        cpu: n_cpu
    }

}

task processResult {
    input {
        File lifted
        File reference
        String? outputName
    }

    command <<<
        python3  /usr/local/compgen/analysis/liftoff/computePeptideLengthFractions.py ~{lifted} ~{reference} > ~{outputName}
    >>>

    runtime {
        docker: "doejgi/compgen@sha256:574ca61c2f96fbfd0ce4623e0d924a3564fa47109fb63fd3125da69bf395eed3"
        cpu: 1
    }

    output {
       File result = "~{outputName}"
    }

}
