### COVID pipeline

A simple pipeline for the assembly, consensus calling and annotation of COVID.

This pipeline is somewhat dependent on folder structures and naming schemes used by the HPC/sequencing core at Icahn School of Medicine.

It also has the primer set hard-coded.

Please ensure you understand how it works and what changes will need to be made if you plan on running on your own data. 
I am happy to help (time permitting) so please reach out.


#### To run this pipeline first create the conda environment for it to run in.

``` cd  <this directory>```

```conda env create --file env.yml```

```conda activate COVID```

you may need to update conda before creation

```conda update conda```

#### Then for pacbio CCS reads run:

```python run_pipeline.py -p <pacbio_ccs_reads.fastq> -o <path/to/output>```

Threads can be specified with ```-t```


#### For Illumina reads run:

```python run_pipeline.py -i <sample_folder>```


The following folder structure should exist
 
```<sample_folder>
└───<reads_15kb_primers>'
│   │   <read_prefix>_1.fastq.gz
│   │   <read_prefix>_2.fastq.gz
│
└───<reads_2kb_primers>
    │   <read_prefix>_1.fastq.gz
    │   <read_prefix>_2.fastq.gz
```

n.b. can be run on one to as many read files as needed, each pair of reads should have it's own folder.

will create **pipeline** folder in <sample_folder> with output.



##### running QC

Finally to run QC do ```python run_QC.py -i <sample_folder> -kb /path/to/kraken_db```

on any <sample_folder> that has run successfully.