import os
import sys
import shutil


if len(sys.argv) != 3:
    sys.stderr.write("""process_run - process an Illumina sars-ncov-2 run.
    USE: looks for folders with a pair of reads in them, iterprets everything before the first hyphen as a sample name
    then copies directories into output folder under sample name.
    USAGE: python process_run.py <illumina_run_dir> <output_dir>""")


for i in os.listdir(sys.argv[1]):
    if os.path.isdir(os.path.join(sys.argv[1], i)):
        sample_name = i.split('-')[0]
        gotr1 = False
        gotr2 = True
        for j in os.listdir(os.path.join(sys.argv[1], i)):
            if j.endswith(".fastq.gz") and j[-15:-13] == "R1":
                gotr1 = True
            if j.endswith(".fastq.gz") and j[-15:-13] == "R2":
                gotr2 = True
        if gotr1 and gotr2 and len(os.listdir(os.path.join(sys.argv[1], i))) == 2:
            sys.stderr.write("Moving %s to sample %s in output directory.\n" % (i, sample_name))
            shutil.copytree(os.path.join(sys.argv[1], i), os.path.join(sys.argv[2], sample_name, i))
        elif not gotr1 or not gotr2:
            sys.stdout.write("Reads not found in folder %s.\n" % i)
        else:
            sys.stderr.write("More than 2 files found in folder %s.\n" % i)
    else:
        sys.stderr.write("%s is not a folder.\n" % i)
