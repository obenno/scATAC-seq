# scATAC-seq

## Pipeline Summary

1. Add barcode to reads (`sinto`)

2. Reads QC and adaptor trimming (`cutadapt`, `fastqc`, `cutqc`)

3. Mapping (`bowtie2`)

4. Generate fragments file (`sinto`)

<details>
<summary>How the fragment file is generated</summary>

From sinto's user guide: https://timoast.github.io/sinto/basic_usage.html

Generating the fragment file involves the following steps in order:

1. Extract cell barcode sequence associated with the fragment.

2. Adjust alignment positions for the 9 bp Tn5 shift by applying +4/-5 to the start and end position of the paired reads.

3. Remove fragments where either read has a MAPQ score less than the specified cutoff.

4. Remove fragments where the fragment size is greater than the specified maximum.

5. Collapse PCR duplicates:

  1. Count the frequency of each fragment for each cell barcode.

  2. Within a cell barcode, collapse fragments that share a start or end coordinate on the same chromosome.

  3. Across all cell barcodes, collapse fragments that share the exact start and end coordinates on the same chromosome.

  4. Assign the fragment to the most abundant cell barcode.

  5. Record the read count for the collapsed fragment.

6. Write fragments to file. Note that fragments are not sorted or compressed.

</details>


5. Library metrics
  - proportion of mitocondrial reads
  - library comlexity (**NRF**, **PBC1**, **PBC2**)
  - fragment length distribution

6. Call peaks as bulk ATAC-seq (`macs2`, `{Signac}`)

7. Cells calling by fragments in peaks

5. Cell QC and filtering (`Signac`)

6. Clustering

7. Post-analysis
