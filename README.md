# trackplot_snakemake

Snakemake pipeline for trackplot to make coverage track plots across a series of intervals

## Running trackplot across a series of intervals

### Ensure all input interval files are compressed with bgzip and are tabix-indexed

Trackplot would do this under the hood anyway. But as we want to immediately parallelise across input intervals, we want to avoid race conditions between different invocations of trackplot. The safest way to do this is to pre-validate the input interval files (GTF file and intervals.txt file). Use `workflow/scripts/check_interval_inputs.py` to do this. Requires samtools installation and python 3.9+:

```bash
python workflow/scripts/check_interval_inputs2.py -g data/gencode.v40.annotation.gtf.gz -i data/intervals.txt -o data/intervals.validated.txt
Running: tabix -p gff data/gencode.v40.annotation.gtf.gz
Error in indexing gff file: Command '['tabix', '-p', 'gff', 'data/gencode.v40.annotation.gtf.gz']' returned non-zero exit status 1.
STDERR: [tabix] the compression of 'data/gencode.v40.annotation.gtf.gz' is not BGZF

Direct indexing failed. The file might not be properly bgzipped. Attempting to decompress and re-bgzip...
Running: sort -k1,1 -k4,4n -o /tmp/tmpx2p9hxyz.sorted /tmp/tmp1rlkdrh7.uncompressed
Created properly bgzipped file at: data/gencode.v40.annotation.gtf.bgz
Running: tabix -p gff data/gencode.v40.annotation.gtf.bgz
Successfully indexed re-bgzipped file: data/gencode.v40.annotation.gtf.bgz
GTF file: data/gencode.v40.annotation.gtf.bgz
Running: sort -k1V -k2n -k3n -o /tmp/tmpdrpwhm2j.sorted data/le_ids.merged.sorted.bed
Created bgzipped file: data/le_ids.merged.sorted.bed.bgz
Running: tabix -p bed data/le_ids.merged.sorted.bed.bgz
Successfully indexed file: data/le_ids.merged.sorted.bed.bgz
Updated intervals file saved to: data/intervals.validated.txt
Intervals file: data/intervals.validated.txt

All input files have now been validated or re-processed
```

Note: only checks file extensions and not the actual file contents themselves

### Dry run:

```bash
snakemake -p --configfile config/config.test.tdp43-apa.yaml --cores 2 --use-singularity --singularity-args="--bind /home/sam" -n
```