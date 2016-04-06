perl /rhome/cjinfeng/software/tools/mirdeep2/bwa_sam_converter.pl -i SRR771499_1.trim3_5.chr2.sam -c -o reads_collapsed.fa -a reads_collapsed_vs_genome.arf
perl ./scripts/tag2Rfam.pl clean.fas ./database/Rfam/Rfam.fasta ./test

perl reformatfasta.pl --fasta hairpin_osa.fa
perl reformatfasta.pl --fasta mature_plant.fa
perl reformatfasta.pl --fasta mature_osa.fa
perl /rhome/cjinfeng/software/bin/getlargeseq.pl --input reads_collapsed.fa --output reads_collapsed.long.fa --length 17
awk '$2>=17' reads_collapsed_vs_genome.arf > reads_collapsed_vs_genome.long.arf

perl /rhome/cjinfeng/software/tools/mirdeep2/miRDeep2.pl reads_collapsed.fa /rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa reads_collapsed_vs_genome.arf
perl /rhome/cjinfeng/software/tools/mirdeep2/miRDeep2.pl reads_collapsed.fa /rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa reads_collapsed_vs_genome.arf ./database/mirbase/mature_osa.fa ./database/mirbase/mature_plant.fa ./database/mirbase/hairpin_osa.fa > log 2> log2 &
perl /rhome/cjinfeng/software/tools/mirdeep2/miRDeep2.pl reads_collapsed.fa /rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa reads_collapsed_vs_genome.arf ./database/mirbase/mature_osa.fa.reform ./database/mirbase/mature_plant.fa.reform ./database/mirbase/hairpin_osa.fa.reform > log 2> log2 &
perl /rhome/cjinfeng/software/tools/mirdeep2/miRDeep2.pl reads_collapsed.long.fa /rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa reads_collapsed_vs_genome.long.arf ./database/mirbase/mature_osa.fa.reform ./database/mirbase/mature_plant.fa.reform ./database/mirbase/hairpin_osa.fa.reform > log 2> log2 &

