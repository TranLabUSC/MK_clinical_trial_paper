source ~/.bashrc
conda activate /project/dtran642_927/SonLe/USC_Source/source/Single_Cell/NextFlow/envs


nextflow run GSEA.nf -c nextflow.GSEA.config --resume
