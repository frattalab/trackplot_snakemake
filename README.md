# trackplot_snakemake

Snakemake pipeline for trackplot to make coverage track plots across a series of intervals

## Installation

Install the execution pipeline using conda/mamba:

```bash
conda env create -f workflow/envs/trackplot_snakemake.yaml
conda activate trackplot_snakemake
```

This minimal environment installs snakemake (>=8.0.0) and samtools (>=1.21.0) to be able to run `bgzip` and `tabix` on input interval files prior to running the pipeline. All snakemake rules themselves are executed inside Singularity containers.

## Running trackplot across a series of intervals

### Ensure all input interval files are compressed with bgzip and are tabix-indexed

Trackplot would do this under the hood anyway. But as we want to immediately parallelise across input intervals, we want to avoid race conditions between different invocations of trackplot. The safest way to do this is to pre-validate the input interval files (GTF file and intervals.txt file). Use `workflow/scripts/check_interval_inputs.py` to do this. Requires samtools installation and python 3.9+ (satisfied by the `trackplot_snakemake` conda environment):

```bash
python workflow/scripts/check_interval_inputs.py -g data/gencode.v40.annotation.gtf.gz -i data/intervals.txt -o data/intervals.validated.txt
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

#### Notes

- Only checks file extensions and not the actual file contents themselves.
- Uses temporary files via the [`tempfile` library](https://docs.python.org/3/library/tempfile.html#tempfile.gettempdir). Use environment variables to modify the temporary directory if necessary.
- Supports BED or GFF/GTF input interval files only (via the tabix presets).

### Configuration

The pipeline requires as input:

- BED files for the 'view' regions (plotting window, passed to `trackplot -e`) and a 'focus' region (highlighted window, passed to `trackplot --focus`). Focus regions are matched to view regions by the 'name' (4th) field of the BED file. The 'focus' region BED file path **must be provided, but it is not required to provide a focus region for any of the input regions**
  - 'view' BED should ideally be BED6 (to include the strand column), the 'focus' bed only needs to be BED4 format (to include the name field)
- density.txt (`trackplot --density`) and intervals.txt (`trackplot --interval`) files as specified by trackplot. This pipeline does not validate the format of the density.txt file in any way. The interval.txt file is checked to ensure the file extensions imply bgzipped and tabix-indexed

Consult [trackplot documentation](https://trackplot.readthedocs.io/en/latest/command/) for instructions on how to generate trackplot-specific input formats.

It's possible to specify multiple datasets per pipeline run/set of input intervals. Simply add an additional dataset_name and density.txt file to the `density_lists` option. Trackplot is run once for each input region and density.txt file. Output PDF files are then concatenated separately **per interval across datasets** (`<region_name>.all_datasets.pdf`) and **per dataset across intervals** (`<dataset_name>.all_regions.pdf`)

### Usage

To perform a dry run on the test data:

```bash
snakemake -p --configfile config/config.test.tdp43-apa.yaml --cores 2 --use-singularity --singularity-args="--bind /home/sam" -n
```