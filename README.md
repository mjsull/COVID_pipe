### COVID pipeline

A simple pipeline for the assembly, consensus calling and annotation of COVID.

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

```python run_pipeline.py -i <illumina_directory> -o <output_folder> -s <samle_name>```


The following folders should exist in ```illumina_directory```

```<sample_name>-1-5p1 <sample_name>-1-5p2 <sample_name>-2p1 <sample_name>-2p2```

or

```<sample_name>-1-5p1-1step <sample_name>-1-5p2-1step <sample_name>-2p1-1step <sample_name>-2p2-1step```

each folder should contain a single pair of reads.