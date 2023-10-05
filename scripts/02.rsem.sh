#!bin/sh
# Author: Hiroshi Kumagai
# Editor: Emma Kumagai
# Transcript quantification with RSEM

PROJECTDIR="/Volumes/HDD14TB/RNAseq_expression/MOTSc_cast_unified"
FILEIDS="${PROJECTDIR}/filename.txt"
INDIR="${PROJECTDIR}/star_rsem"
OUTDIR="${INDIR}/rsem_results"
INDEXDIR="/Volumes/HDD14TB/RNAseq_expression/Tools/Ensembl_RSEM-1.3.3_NUC_index_mus/Ensembl_RSEM_Mus_index"

cat ${FILEIDS} | while read line; do
    rsem-calculate-expression \
        --alignments \
        --strandedness reverse \
        --estimate-rspd \
        --no-bam-output \
        -p 12 \
        ${INDIR}/${line}.Aligned.toTranscriptome.out.bam \
        ${INDEXDIR} \
        ${OUTDIR}/${line}
done
