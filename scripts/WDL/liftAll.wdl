version 1.0

import "liftoff.wdl" as liftoff

workflow liftAll {

    input {
        Array[Pair[File, File]] genomeAnnotationPairs
        Int liftoff_n_cpu = 16
    }

    Array[Int] indices = range(length(genomeAnnotationPairs))
    scatter (i in indices) {
        scatter (j in indices) {
            if (i != j) {
                String refGenome = basename(genomeAnnotationPairs[j].left)
                call liftoff.liftoff {
                    input:  targetFasta=genomeAnnotationPairs[i].left,
                            referenceFasta=genomeAnnotationPairs[j].left,
                            referenceGff=genomeAnnotationPairs[j].right,
                            outputName="~{refGenome}_to_~{targetGenome}.gff3",
                            n_cpu=liftoff_n_cpu
                }
            }
        }

        String targetGenome = basename(genomeAnnotationPairs[i].left)
        Array[File] targetResults = select_all(liftoff.result)
        call catResults {
            input:
                files=targetResults,
                outName="liftoffs_to_~{targetGenome}.gff3"
        }
        Pair[String, File] targetMap = (targetGenome, catResults.result)
    }

    output {
        Array[Pair[String, File]] resultsByTarget  = targetMap
    }
}

task catResults {
    input {
        Array[File] files
        String      outName
    }

    command <<<
        set -e
        echo "##gff-version 3" > ~{outName}
        
        for file in ~{sep=' ' files}; do
            grep -v "^#" $file >> ~{outName}
        done
    >>>

    output {
        File result = outName
    }

    runtime {
        docker: "debian@sha256:98d3b4b0cee264301eb1354e0b549323af2d0633e1c43375d0b25c01826b6790"
        cpu: 1
        memory: "4G"
    }
}
