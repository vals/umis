rm test*.fq
rm test*.sam

umis fastqtransform \
examples/MARS-Seq/transform.json \
examples/MARS-Seq/SRP035326.fastq \
> test01.fq

umis fastqtransform \
examples/MARS-Seq/transform.json \
examples/MARS-Seq/SRP035326_5.fastq \
> test02.fq

umis fastqtransform \
examples/CEL-Seq/transform.json \
examples/CEL-Seq/SRP036633_1.fastq \
examples/CEL-Seq/SRP036633_2.fastq \
> test03.fq

umis fastqtransform \
examples/DropSeq/transform.json \
examples/DropSeq/SRR1873278_1.fastq \
examples/DropSeq/SRR1873278_2.fastq \
> test04.fq

umis fastqtransform \
examples/inDrop/transform.json \
examples/inDrop/SRR1784317_1.fastq \
examples/inDrop/SRR1784317_2.fastq \
> test05.fq

umis fastqtransform \
--demuxed_cb ATTAGAC \
examples/STRT-Seq/SRP022764_transform.json \
examples/STRT-Seq/SRP022764_ESCell_1_ATTAGAC_single.fastq \
> test06.fq

umis fastqtransform \
--demuxed_cb A01 \
examples/STRT-Seq/SRP045452_transform.json \
examples/STRT-Seq/SRP045452_1772058148_A01.fastq \
> test07.fq

umis fastqtransform \
--demuxed_cb CACTGT \
examples/BATSeq/transform.json \
examples/BATseq/SRR1558183_1.fastq \
examples/BATseq/SRR1558183_2.fastq \
> test08.fq

umis fastqtransform \
examples/CEL-Seq/transform.json \
examples/CEL-Seq/SRP036633_1.fastq.gz \
examples/CEL-Seq/SRP036633_2.fastq.gz \
> test09.fq

umis fastqtransform \
examples/CEL-Seq/transform.json \
examples/CEL-Seq/SRP048838_1.fastq \
examples/CEL-Seq/SRP048838_2.fastq \
> test10.fq

umis fastqtransform \
--dual_index \
examples/STRT-Seq/dual_index_transform.json \
examples/STRT-Seq/dualindex_example_1.fastq \
examples/STRT-Seq/dualindex_example_2.fastq \
> test11.fq

umis fastqtransform \
--keep_fastq_tags \
--umi_only \
--fastq1out test12_1.fq \
--fastq2out test12_2.fq \
examples/paired-with-umi-read/transform.json \
examples/paired-with-umi-read/fq_1.fq \
examples/paired-with-umi-read/fq_2.fq \
examples/paired-with-umi-read/umi.fq 

umis bamtag \
examples/bamtag/bamtag.sam \
> test_bamtag.sam
