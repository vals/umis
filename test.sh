rm -r tests/results
mkdir -p tests/results

umis fastqtransform \
examples/MARS-Seq/transform.json \
examples/MARS-Seq/SRP035326.fastq \
> tests/results/test01.fq

umis fastqtransform \
examples/MARS-Seq/transform.json \
examples/MARS-Seq/SRP035326_5.fastq \
> tests/results/test02.fq

umis fastqtransform \
examples/CEL-Seq/transform.json \
examples/CEL-Seq/SRP036633_1.fastq \
examples/CEL-Seq/SRP036633_2.fastq \
> tests/results/test03.fq

umis fastqtransform \
examples/DropSeq/transform.json \
examples/DropSeq/SRR1873278_1.fastq \
examples/DropSeq/SRR1873278_2.fastq \
> tests/results/test04.fq

umis fastqtransform \
examples/inDrop/transform.json \
examples/inDrop/SRR1784317_1.fastq \
examples/inDrop/SRR1784317_2.fastq \
> tests/results/test05.fq

umis fastqtransform \
--demuxed_cb ATTAGAC \
examples/STRT-Seq/SRP022764_transform.json \
examples/STRT-Seq/SRP022764_ESCell_1_ATTAGAC_single.fastq \
> tests/results/test06.fq

umis fastqtransform \
--demuxed_cb A01 \
examples/STRT-Seq/SRP045452_transform.json \
examples/STRT-Seq/SRP045452_1772058148_A01.fastq \
> tests/results/test07.fq

umis fastqtransform \
--demuxed_cb CACTGT \
examples/BATSeq/transform.json \
examples/BATseq/SRR1558183_1.fastq \
examples/BATseq/SRR1558183_2.fastq \
> tests/results/test08.fq

umis fastqtransform \
examples/CEL-Seq/transform.json \
examples/CEL-Seq/SRP036633_1.fastq.gz \
examples/CEL-Seq/SRP036633_2.fastq.gz \
> tests/results/test09.fq

umis fastqtransform \
examples/CEL-Seq/transform.json \
examples/CEL-Seq/SRP048838_1.fastq \
examples/CEL-Seq/SRP048838_2.fastq \
> tests/results/test10.fq

umis fastqtransform \
examples/STRT-Seq/dual_index_transform.json \
examples/STRT-Seq/dualindex_example_1.fastq \
examples/STRT-Seq/dualindex_example_2.fastq \
> tests/results/test11.fq

umis fastqtransform \
--keep_fastq_tags \
--fastq1out tests/results/test12_1.fq \
--fastq2out tests/results/test12_2.fq \
examples/paired-with-umi-read/transform.json \
examples/paired-with-umi-read/fq_1.fq \
examples/paired-with-umi-read/fq_2.fq \
examples/paired-with-umi-read/umi.fq

umis fastqtransform \
examples/SCRB-Seq/transform.json \
examples/SCRB-Seq/scrbseq_R1.fastq \
examples/SCRB-Seq/scrbseq_R2.fastq \
> tests/results/test13.fq

umis fastqtransform \
--separate_cb \
examples/Klein-inDrop/transform.json \
examples/Klein-inDrop/klein-v3_R1.fq \
examples/Klein-inDrop/klein-v3_R2.fq \
examples/Klein-inDrop/klein-v3_R3.fq \
examples/Klein-inDrop/klein-v3_R4.fq \
> tests/results/test14.fq

umis demultiplex_samples --nedit 1 \
--barcodes examples/Klein-inDrop/sample-index.txt \
--out_dir tests/results \
examples/Klein-inDrop/test14.fq

umis bamtag \
examples/bamtag/bamtag.sam \
> tests/results/test_bamtag.sam

# test streaming bamtag
umis bamtag - < \
examples/bamtag/bamtag.bam \
> tests/results/test_streaming_bamtag.sam

# test conflicting tag/qname annotations
umis bamtag - < examples/bamtag/bamtag-xstag.bam > tests/results/test-bamtag-xstag.sam

umis cb_histogram \
examples/Klein-inDrop/test14.fq \
> tests/results/test15-cb-histogram.txt

umis cb_histogram \
--umi_histogram tests/results/test15-mb-histogram.txt \
examples/Klein-inDrop/test14.fq \
> tests/results/test15-cb-histogram.txt

umis umi_histogram \
examples/Klein-inDrop/test14.fq \
     > tests/results/test16-umi-histogram.txt

umis tagcount \
     examples/tagcount/tagcount.sam \
     tests/results/test17-tagcount.txt

# only display diff output if there are differences
if [[ $(diff -rq tests/results tests/correct) ]]; then
  diff -rq tests/results tests/correct
  exit 1
else
  echo "Tests passed."
fi
