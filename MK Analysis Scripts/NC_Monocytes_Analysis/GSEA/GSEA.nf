#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process GSEA_analysis {
    clusterOptions " --nodes=1 --ntasks=1 --mem-per-cpu=64GB --cpus-per-task=1 --time=3:00:00 --partition=epyc-64"
    conda 'java'

    input:
    path rnk
    each geneset
    val geneset_folder
    output:


    script:
    """
    #!/bin/bash

    # Get the real path of the rnk file
    real_path=\$(readlink -f ${rnk})

    echo rnk_file: \${real_path}

    # Create output directory one level up from .rnk file
    PA_dir=\$(dirname \${real_path})/../PathwayAnalysis
    [ ! -d "\${PA_dir}" ] && mkdir \${PA_dir}

    IFS="=" read -r geneset_name geneset_file <<< "$geneset"

    outdir=\$PA_dir/\${geneset_name}
    [ ! -d "\${outdir}" ] && mkdir \${outdir}
    echo outdir:\${outdir}
    echo geneset:\${geneset_name}
    echo geneset_folder: ${geneset_folder}

    source ~/.bashrc

    \${SOURCE}/GSEA/gsea-cli.sh GSEAPreranked \
    -gmx ${geneset_folder}/\${geneset_file} \
    -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 \
    -rnk ${rnk}  \
    -scoring_scheme classic -rpt_label my_analysis \
    -create_svgs true  \
    -include_only_symbols true \
    -make_sets true -plot_top_x 300 \
    -rnd_seed timestamp -set_max 1000 \
    -set_min 5 -zip_report false \
    -out \${outdir}
    """
}

workflow {
    // Find .rnk files and create a channel
    rnk_files = Channel.fromPath("${params.indir}/**/*.rnk", type: 'file')
    // rnk_files = Channel.fromPath("${params.indir}/*.rnk", type: 'file')
    // Print each file in the rnk_files channel
    rnk_files.view { "rnk file: $it" }
    GSEA_analysis(rnk_files, params.genesets, params.geneset_folder)
}
