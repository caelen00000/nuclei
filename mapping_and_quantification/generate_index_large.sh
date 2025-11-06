#!/bin/bash

STAR \
--runThreadN 36 \
--runMode genomeGenerate \
--genomeDir genome_dir_large \
--genomeFastaFiles "/mnt/z/Caelen/Resources/Flybase/dmel-all-chromosome-r6.62.fasta" \
--sjdbGTFfile "/mnt/z/Caelen/Resources/Flybase/dmel-all-r6.62.gtf" \
--genomeSAindexNbases 12 \
--sjdbOverhang 119

# --genomeSAindexNbases 12 is just an adjustment to the algorithm due to the relatively small genome of drosophila
# --sjdbOverhang 119 this is the maximum read length - 1