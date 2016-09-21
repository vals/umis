# umis


**umis** provides tools for estimating expression in RNA-Seq data which performs
sequencing of end tags of trancsript, and incorporate molecular tags to
correct for amplification bias.

There are four steps in this process.

 1. Formatting reads
 2. Filtering noisy cellular barcodes
 3. Pseodomapping to cDNAs
 4. Counting molecular identifiers

## 1. Formatting reads

We want to strip out all non-biological segments of the sequenced reads for
the sake of mapping. While also keeping this information for later use. We
consider non-biological information such as Cellular Barcode and Molecular
Barcode. To later be able to extract the optional CB and the MB these are put
in the read header, with the followign format.

    @HWI-ST808:130:H0B8YADXX:1:1101:2088:2222:CELL_GGTCCA:UMI_CCCT
    AGGAAGATGGAGGAGAGAAGGCGGTGAAAGAGACCTGTAAAAAGCCACCGN
    +
    @@@DDBD>=AFCF+<CAFHDECII:DGGGHGIGGIIIEHGIIIGIIDHII#

The command `umis fastqtransform` is for transforming a (pair of) read(s) to
this format based on a _transform file_. The transform file is a json file
which has a Python flavored regular expression for each read, made to extract
the necessary components of the reads.

## 2. Filtering noisy cellular barcodes
Not all cellular barcodes identified in the transformation will be real. Some
will be low abundance barcodes that do not represent an actual cell. Others
will be barcodes that don't come from a set of known barcodes. The `umi cb_filter`
command can be used to filter a transformed FASTQ file, dropping unknown
barcodes. The `--nedit` option can be supplied to correct barcodes `--nedit`
distance away from known barcodes. After barcode filtering,
the `umi cb_histogram` command will generate a file of counts for
each cellular barcode. This file can be used to find a count cutoff for barcodes
that are high abundance for downstream quantitation.

## 3. Pseudomapping to cDNAs

This is done by pseduoaligners, either Kallisto or RapMap. The SAM file output
from these tools need to be saved.

## 4. Counting molecular identifiers

The final step is to infer which cDNA was the origin of the tag a UMI was
attached to. We use the pseudoalignments to the cDNAs, and consider a tag
assigned to a cDNA as a partial _evidence_ for a (cDNA, UMI) pairing. For
actual counting, we only count unique UMIs for (gene, UMI) pairings with
sufficient evidence.

## kallisto
The quantitation used in `umis` handles reads that could come from multiple
transcripts by assigning a fractional count to each transcript and then
filtering for a minimum count at the end. Many single-cell analyses use
something similar to this type of counting, but it has drawbacks
(see
[this paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0970-8)).
For more principled UMI quantification,
see [Kallisto](https://github.com/pachterlab/kallisto). kallisto needs the files
in a certain format: each cellular barcode has its own FASTQ file and a file
that lists the UMI for each read. The `umis kallisto` command can reformat your
fastq files to that format.
