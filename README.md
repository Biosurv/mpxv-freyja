# mpxv-freyja

**This pipeline is still in active development**

Simple reference-based assembly of Oxford Nanopore data (amplicon schemes, IDT xGen panels or
primer-free capture methods such as Twist) with medaka, built for wastewater MPXV method testing. It reports read, mapping
and coverage statistics and, optionally, Freyja clade abundances.

## Requirements

* Linux with Nextflow ≥ 26.04 (`conda install -c bioconda nextflow`)
* Docker, Singularity/Apptainer or conda

## Run

```bash
nextflow run mpxv-freyja -profile docker \\
    --reads /data/fastq\_pass \\
    --reference reference.fasta \\
    --primers primers.bed \\
    --freyja\_barcodes barcodes.csv \\
    --outdir results
```

## Inputs and parameters

|Parameter|Default|Description|
|-|-|-|
|`--reads`|required|Folder with one sub-folder of FASTQ files per sample (e.g. MinKNOW `fastq\_pass`); the sub-folder name is the sample name. FASTQ files placed directly in the folder are treated as one sample each. `unclassified` is skipped|
|`--reference`|required|Reference FASTA (the one the primer scheme was designed on, when using primers)|
|`--primers`|–|Primer BED (`<name>\_<n>\_LEFT/RIGHT`) or IDT xGen masterfile; the format is detected automatically. Leave out for primer-free methods (e.g. Twist): no primer trimming is done and the target/amplicon columns are left blank|
|`--freyja\_barcodes`|–|Freyja barcode CSV|
|`--freyja\_lineages`|–|Optional lineage hierarchy YAML for the barcodes|
|`--min\_depth`|`20`|Depth threshold for coverage statistics; positions below it are N in the consensus|
|`--min\_read\_length`|`0`|Minimum read length (`0` = no limit)|
|`--max\_read\_length`|`0`|Maximum read length (`0` = no limit)|
|`--medaka\_model`|`auto`|medaka model; `auto` reads the basecaller model from the FASTQ headers|
|`--outdir`|`results`|Output folder|

## Outputs

|File|Contents|
|-|-|
|`summary.csv`|One row per sample. Reads: total, after length filter, mapped, after primer trim. Depth threshold. Genome and target length, mean depth and % covered at ≥ `min\_depth`. Amplicons below threshold. Consensus length and Ns, SNPs and indels, medaka model. With Freyja: status, clade summary, top lineage and % per clade|
|`amplicon\_coverage.csv`|Mean depth and % of each amplicon insert covered at ≥ `min\_depth`, per sample (only with `--primers`)|
|`freyja\_lineages.csv`|Freyja lineage abundances per sample, with the clade taken from the lineage name (`MPX-IIb.C.1` → `IIb`)|
|`consensus/<sample>.consensus.fasta`|Consensus sequence|

## Freyja barcodes

* The barcode positions must use the coordinates of `--reference`. The pipeline checks this, and
if they don't match it reports `skipped`.

