#!/bin/sh
# 1. Align the reads with STAR (https://github.com/alexdobin/STAR)
# 2. Quantify the expression with RSEM (https://github.com/deweylab/RSEM)

SEQLIBS=(SRR6763185 SRR6763193 SRR6763201 SRR6763209 SRR6763217 SRR6763225 SRR6763233 SRR6763241 SRR6763249 SRR6763257 SRR6763265 SRR6763273 SRR6763281 SRR6763289 SRR6763297 SRR6763305 SRR6763313 SRR6763321 SRR6763329)

# STAR
for seqlib in ${SEQLIBS[@]}; do
	STAR \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode TranscriptomeSAM \
		--runThreadN 4 \
		--outSAMattributes All \
		--readFilesCommand gunzip -c \
		--genomeDir /Volumes/HDD14TB/RNAseq_expression/Tools/Ensembl_STAR_Homo_index \
		--readFilesIn ${seqlib}_1.fastq.gz ${seqlib}_2.fastq.gz \
		--outFileNamePrefix /Volumes/HDD14TB/RNAseq_expression/Dataset/Sarcopenia_Chinese_new/star_rsem/${seqlib}.
done


# RSEM
SEQLIBS=(SRR6763185 SRR6763193 SRR6763201 SRR6763209 SRR6763217 SRR6763225 SRR6763233 SRR6763241 SRR6763249 SRR6763257 SRR6763265 SRR6763273 SRR6763281 SRR6763289 SRR6763297 SRR6763305 SRR6763313 SRR6763321 SRR6763329)
for seqlib in ${SEQLIBS[@]}; do
	rsem-calculate-expression \
		--paired-end --alignments \
		--strandedness reverse \
		--estimate-rspd \
		--no-bam-output \
		-p 12 \
		${seqlib}.Aligned.toTranscriptome.out.bam \
		/Volumes/HDD14TB/RNAseq_expression/Tools/Ensembl_RSEM-1.3.3_NUC_index_Homo/Ensembl_RSEM-1.3.3_NUC_index_Homo \
		/Volumes/HDD14TB/RNAseq_expression/Dataset/Sarcopenia_Chinese_new/star_rsem/rsem_results/${seqlib}
done

