STAR --runThreadN 16    \
--readFilesCommand zcat  \
--quantMode TranscriptomeSAM GeneCounts   \
--genomeDir ~/data/genome/eb2020/STAR_db/   \
--readFilesIn ./B6R3_FRAS210078503-1r_1.clean.fq ./B6R3_FRAS210078503-1r_2.clean.fq   \
--outSAMattrRGline ID:B6R3 SM:B6R3 PL:ILLUMINA    \
--outSAMtype BAM SortedByCoordinate   \
--outFileNamePrefix B6R3   \
--outReadsUnmapped Fastx  