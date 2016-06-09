# umis

**Note: This tool works on heuristic counting. For more principled UMI quantification, see [Kallisto](https://github.com/pachterlab/kallisto).**
Some scripts of `umis` might still be useful for pre-processing and investigating UMI data. In particular the regex based UMI extraction.

**umis** provides tools for estimating expression in RNA-Seq data which performs
sequencing of end tags of trancsript, and incorporate molecular tags to
correct for amplification bias.

There are three steps in this process.

 1. Formatting reads
 2. Pseodomapping to cDNAs
 3. Counting molecular identifiers

## 1. Formatting reads

We want to strip out all non-biological segments of  the sequenced reads for
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

## 2. Pseodomapping to cDNAs

This is done by pseduoaligners, either Kallisto or RapMap. The SAM file output
from these tools need to be saved.

## 3. Counting molecular identifiers

The final step is to infer which cDNA was the origin of the tag a UMI was
attached to. We use the pseudoalignments to the cDNAs, and consider a tag
assigned to a cDNA as a partial _evidence_ for a (cDNA, UMI) pairing. For
actual counting, we only count unique UMIs for (gene, UMI) pairings with
sufficient evidence.
