### COVID pipeline

A simple pipeline for the assembly, consensus calling and annotation of COVID.

To run this pipeline first create the conda environment for it to run in.

``` cd  <this directory>```
```conda env create --file env.yml```
```conda activate COVID```

you may need to update conda before creation

```conda update conda```

Then for pacbio CCS reads run:

```python run_pipeline.py -p <pacbio_ccs_reads.fastq> -o <path/to/output>```

Threads can be specified with ```-t```
