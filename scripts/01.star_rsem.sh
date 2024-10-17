# File names
SEQLIBS=(BJ_100001 BJ_100002 BJ_100003 BJ_100004 BJ_100005 BJ_100006 BJ_100007 BJ_100008 BJ_100009 BJ_100011 BJ_100012 BJ_100013 BJ_100014 BJ_100015 BJ_100016 BJ_100017 BJ_100018 BJ_100019 BJ_100020 BJ_100021 BJ_100022 BJ_100023 BJ_100024 BJ_100025 BJ_100026 BJ_100027 BJ_100028 BJ_100029 BJ_100030 BJ_100031 BJ_100032 BJ_100033 BJ_100034 BJ_100035 BJ_100036 BJ_100037 BJ_100038 BJ_100039 BJ_100040 BJ_100041 BJ_100042 BJ_100043 BJ_100044 BJ_100045 BJ_100046 BJ_100047 BJ_100048 BJ_100049)

# Runing STAR
for seqlib in ${SEQLIBS[@]}; do
	STAR \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode TranscriptomeSAM \
		--runThreadN 4 \
		--outSAMattributes All \
		--readFilesCommand gunzip -c \
		--genomeDir /Volumes/HDD14TB/RNAseq_expression/Tools/Ensembl_STAR_Homo_index \
		--readFilesIn ${seqlib}_1.fastq.gz ${seqlib}_2.fastq.gz \
		--outFileNamePrefix /Volumes/HDD14TB/RNAseq_expression/Dataset/biopsy/${seqlib}.
done
