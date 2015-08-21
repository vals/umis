rm test*.fq

umis fastqtransform \
examples/MARS-Seq/transform.json \
examples/MARS-Seq/SRP035326.fastq \
> test1.fq

umis fastqtransform \
examples/MARS-Seq/transform.json \
examples/MARS-Seq/SRP035326_5.fastq \
> test2.fq

umis fastqtransform \
examples/CEL-Seq/transform.json \
examples/CEL-Seq/SRP036633_1.fastq \
examples/CEL-Seq/SRP036633_2.fastq \
> test3.fq

umis fastqtransform \
examples/DropSeq/transform.json \
examples/DropSeq/SRR1873278_1.fastq \
examples/DropSeq/SRR1873278_2.fastq \
> test4.fq

umis fastqtransform \
examples/inDrop/transform.json \
examples/inDrop/SRR1784317_1.fastq \
examples/inDrop/SRR1784317_2.fastq \
> test5.fq

umis fastqtransform \
--demuxed_cb ATTAGAC \
examples/STRT-Seq/SRP022764_transform.json \
examples/STRT-Seq/SRP022764_ESCell_1_ATTAGAC_single.fastq \
> test6.fq

umis fastqtransform \
--demuxed_cb A01 \
examples/STRT-Seq/SRP045452_transform.json \
examples/STRT-Seq/SRP045452_1772058148_A01.fastq \
> test7.fq

umis fastqtransform \
--demuxed_cb CACTGT \
examples/BATSeq/transform.json \
examples/BATseq/SRR1558183_1.fastq \
examples/BATseq/SRR1558183_2.fastq \
> test8.fq

umis fastqtransform \
examples/CEL-Seq/transform.json \
examples/CEL-Seq/SRP036633_1.fastq.gz \
examples/CEL-Seq/SRP036633_2.fastq.gz \
> test9.fq
