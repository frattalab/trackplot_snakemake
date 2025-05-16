# trackplot_snakemake
Snakemake pipeline for trackplot to make coverage track plots across a series of intervals

## Running trackplot across a series of intervals

dry run:

```bash
snakemake -p --configfile config/config.test.tdp43-apa.yaml --cores 2 --use-singularity --singularity-args="--bind /home/sam" -n
```