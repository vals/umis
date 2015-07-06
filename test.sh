python umicount.py fastqtransform \
--fastq1 examples/MARS-Seq/SRP035326.fastq \
--transform examples/MARS-Seq/transform.json \
--outfastq test1.fq

python umicount.py fastqtransform \
--fastq1 examples/CEL-Seq/SRP036633_1.fastq \
--fastq2 examples/CEL-Seq/SRP036633_2.fastq \
--transform examples/CEL-Seq/transform.json \
--outfastq test2.fq

python umicount.py fastqtransform \
--fastq1 examples/DropSeq/SRR1873278_1.fastq \
--fastq2 examples/DropSeq/SRR1873278_2.fastq \
--transform examples/DropSeq/transform.json \
--outfastq test3.fq

python umicount.py fastqtransform \
--fastq1 examples/inDrop/SRR1784317_1.fastq \
--fastq2 examples/inDrop/SRR1784317_2.fastq \
--transform examples/inDrop/transform.json \
--outfastq test4.fq
