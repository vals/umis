## 1.0.4 (in progress)
- Enable cb_histogram to be used on samples without UMIs.
- Enable filtering of cells during `demultiplex_cells`.

## 1.0.3 
- Python 3 support

## 1.0.2 
- Add `demultiplex_cells` subcommand to break a transformed FASTQ file into separate FASTQ files by cell.
- Future proofing for changes to panda's `to_csv` function.

## 1.0.1 
- Add support for click 7.0.

## 1.0.0
- Fix for min-length filtering with paired samples. Previously required only one read to be longer, fix requires both.
- Fix tests for fastqtagcount to use indexed BAM files.
- Support gzipped cellular barcode files.
- Support 10x V2 barcoding scheme. Thanks to @tomasgomes for the fix.
- Re-enable streaming for cellular barcode filtering.
- Add `--umi_matrix` option to fasttagcount. This outputs a non-umi-deduped matrix of counts, useful for QC.
- Support gzipped files for `sb_filter`, `mb_filter` and `add_uid`.

## 0.8.0
- Fix `fasttagcount` off-by-one issue.
- Add `version` subcommand.
- Fix missing pandas import in `sparse` subcommand.

## 0.7.0
- Fix for kallisto output failing due to defaultdict not being imported. Thanks to @andreas-wilm for the fix.
- Added `tagcount` option `--parse_tags` to use BAM tags rather than parsing read names (`UM` for UMI, `CR` for cell barcode)
- Added `tagcount` option `--gene_tags` to use BAM tags to get ID of mapping gene (`GX` tag).
- Fix tagcount with `--genemap` option not including a column name for the index.
- Add `sparse` subcommand to turn a matrix into a sparse matrix.
- Add `fasttagcount` subcommand. This assumes the input BAM/SAM file is coordinate sorted. Reduces memory usage by over
  100x and runtime by 30-40% for deep samples.
- Warn, don't fail if transcripts are missing from the genemap. 

## 0.6.0
- Fix skipping first piece of evidence when tagcounting.
- Add test for tagcount.
- Output full sorted transcript table from tagcount rather than only the observed transcripts.
- Add `--sparse` option to output tagcount matrices in MatrixMarket format.
- Allow cb_histogram subcommand to take gzipped files.
- Allow cb_filter subcommand to take gzipped files.
- Add support for triple-cellular barcodes.
- Add example for Illumina SureCell (https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/surecell-wta-ddseq.html)

## 0.5.0

- Fix automatic format detection in cb_histogram.
- Add tests for cb_histogram.
- Re-enable streaming bamtagging. Thanks to @chapmanb for the suggestion.
- Add subset_bamfile to subset a BAM file to keep alignments with a given set of cellular barcodes.
- Speed improvements for reading gzipped FASTQ files.
- Memory usage improvements for tagcount.

## 0.4.0

- Fix for handling unicode, thanks to @chapmanb and @sowmyaiyer
- Adds support for adding BAM tags to aligned fastqtransformed files. Thanks to @chapmanb.
- Adds support for UMI-only fastqtransformation.
- Adds support for paired-end target sequences.
- Adds support for detecting sample barcodes via the SB tag in the regex.
- Adds support for sample-based demultiplexing with error correction.

## 0.3.0

- Now supports transforming 3-file input, as from the Linnarsson lab STRT-Seq data
- New kallisto subcommand formats read files for input to kallisto's UMI mode
- Fix gzip based fastq reading on Python 3.5
- Including preliminary subcommand for guessing cell cutoff from cb_histogram

## 0.2.2

- Added MANIFEST file which broke pip installation
