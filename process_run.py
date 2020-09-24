import os
import sys
import shutil
import argparse




def process_reads(args):

    for i in os.listdir(args.input_folder):
        if os.path.isdir(os.path.join(args.input_folder, i)):
            sample_name = i.rstrip('/')
            gotr1 = False
            gotr2 = True
            for j in os.listdir(os.path.join(args.input_folder, i)):
                if j.endswith(args.read1_suffix):
                    gotr1 = True
                if j.endswith(args.read2_suffix):
                    gotr2 = True
            if gotr1 and gotr2 and len(os.listdir(os.path.join(args.input_folder, i))) == 2:
                sys.stderr.write("Moving %s to sample %s in output directory.\n" % (i, sample_name))
                shutil.copytree(os.path.join(args.input_folder, i), os.path.join(args.output_folder, sample_name, os.path.basename(args.input_folder)))
            elif not gotr1 or not gotr2:
                sys.stdout.write("Reads not found in folder %s.\n" % i)
            else:
                sys.stderr.write("More than 2 files found in folder %s.\n" % i)
        else:
            sys.stderr.write("%s is not a folder.\n" % i)


__version__ = "0.1.1"
parser = argparse.ArgumentParser(prog='COVID pipeline', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description="""process_run - process an Illumina sars-ncov-2 run.
    USE: looks for folders with a pair of reads in them, interpre   ts everything before the first underscore as a sample name
    then copies directories into output folder under sample name.
    USAGE: python process_run.py <illumina_run_dir> <output_dir>""")



parser.add_argument('-i', '--input_folder', action='store', help='place where illumina reads are located')
parser.add_argument('-o', '--output_folder', action='store', default="temp", help='output directory')
parser.add_argument('-r1', '--read1_suffix', action='store', default="R1_001.fastq.gz", help='suffix for finding read 1')
parser.add_argument('-r2', '--read2_suffix', action='store', default="R2_001.fastq.gz", help='suffix for finding read 2')
parser.add_argument('-v', '--version', action='store_true', help="print version and exit")


args = parser.parse_args()

if args.version:
    sys.stdout.write('Version %s' % __version__)
    sys.exit()

process_reads(args)