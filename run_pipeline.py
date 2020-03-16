#!/usr/bin/env python3

import argparse
import subprocess
import sys, os


repo_dir = sys.path[0]
virus_length = 29760

def run_illumina(args):
    subprocess.Popen("cutadapt -j %s -g file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -a "
                     "file:%s/db/SARS-CoV-2_primers_3prime_NI.fa -o %s/reads.1.fq.gz %s > %s/cutadapt.1.log"
                     % (args.threads, repo_dir, repo_dir, args.working_dir, args.read_1, args.working_dir), shell=True).wait()
    subprocess.Popen("cutadapt -j %s -g file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -a "
                     "file:%s/db/SARS-CoV-2_primers_3prime_NI.fa -o %s/reads.2.fq.gz %s > %s/cutadapt.2.log"
                     % (args.threads, repo_dir, repo_dir, args.working_dir, args.read_2, args.working_dir), shell=True).wait()
    with open(args.working_dir + '/cutadapt.1.log') as f:
        for line in f:
            if line.startswith("Total written (filtered):"):
                bp = int(line.split()[3].replace(',', ''))
                break
    with open(args.working_dir + '/cutadapt.2.log') as f:
        for line in f:
            if line.startswith("Total written (filtered):"):
                bp += int(line.split()[3].replace(',', ''))
                break
    if bp / virus_length > args.coverage_pilon:
        downsample = args.coverage_pilon / (bp / virus_length)
        subprocess.Popen("seqtk sample %s/reads.1.fq.gz %f | gzip > %s/reads.pilon.1.fq.gz"
                         % (args.working_dir, downsample, args.working_dir), shell=True).wait()
        subprocess.Popen("seqtk sample %s/reads.1.fq.gz %f | gzip > %s/reads.pilon.1.fq.gz"
                         % (args.working_dir, downsample, args.working_dir), shell=True).wait()
        pilon_read_1 = "%s/reads.pilon.1.fq.gz" % args.working_dir
        pilon_read_2 = "%s/reads.pilon.1.fq.gz" % args.working_dir
    else:
        pilon_read_1 = "%s/reads.1.fq.gz" % args.working_dir
        pilon_read_2 = "%s/reads.2.fq.gz" % args.working_dir
    subprocess.Popen("minimap2 -t %s -ax sr %s/db/COVID.fa %s %s | samtools view -b | samtools sort -@ %s -o %s/ref.bam -"
                     " && samtools index %s/ref.bam"
                     % (args.threads, repo_dir, pilon_read_1, pilon_read_2, args.threads, args.working_dir, args.working_dir), shell=True).wait()
    subprocess.Popen("pilon --fix bases --threads %s --mindepth 20 --genome %s/db/COVID.fa --frags %s/ref.bam --tracks --output %s/pilon"
                     % (args.threads, repo_dir, args.working_dir, args.working_dir), shell=True).wait()
    subprocess.Popen("shovill --outdir %s/shovill --R1 %s --R2 %s --gsize %d --cpus %s"
                     % (args.working_dir, pilon_read_1, pilon_read_2, virus_length, args.threads), shell=True).wait()
    subprocess.Popen("prokka  --cpus %s --outdir %s/prokka --prefix sequence --kingdom Viruses --proteins %s/db/COVID.gbk  %s/pilon.fasta "
                     % (args.threads, args.working_dir, repo_dir, args.working_dir), shell=True).wait()



def run_ccs(args):
    subprocess.Popen("cutadapt -j %s -g file:%s/db/SARS-CoV-2_primers_5prime_anchored.fa -a "
                     "file:%s/db/SARS-CoV-2_primers_3prime_anchored.fa -o %s/reads.1.fq.gz %s > %s/cutadapt.1.log"
                     % (args.threads, repo_dir, repo_dir, args.working_dir, args.ccs_reads, args.working_dir), shell=True).wait()
    with open(args.working_dir + '/cutadapt.1.log') as f:
        for line in f:
            if line.startswith("Total written (filtered):"):
                bp = int(line.split()[3].replace(',', ''))
                break
    if bp / virus_length > args.coverage_pilon:
        downsample = args.coverage_pilon / (bp / virus_length)
        subprocess.Popen("seqtk sample %s/reads.1.fq.gz %f | gzip > %s/reads.pilon.1.fq.gz"
                         % (args.working_dir, downsample, args.working_dir), shell=True).wait()
        pilon_read_1 = "%s/reads.pilon.1.fq.gz" % args.working_dir
    else:
        pilon_read_1 = "%s/reads.1.fq.gz" % args.working_dir
    subprocess.Popen("minimap2 -t %s -ax map-pb %s/db/COVID.fa %s | samtools view -b | samtools sort -@ %s -o %s/ref.bam -"
                     " && samtools index %s/ref.bam"
                     % (args.threads, repo_dir, pilon_read_1, args.threads, args.working_dir, args.working_dir), shell=True).wait()
    subprocess.Popen("pilon --fix bases --threads %s --mindepth 20 --genome %s/db/COVID.fa --unpaired %s/ref.bam --tracks --output %s/pilon"
                     % (args.threads, repo_dir, args.working_dir, args.working_dir), shell=True).wait()
    with open(args.working_dir + '/pilon.fasta') as f:
        seq = ''
        for line in f:
            if not line.startswith('>'):
                seq += line.rstrip()
    seq = list(seq)
    with open(args.working_dir + '/pilonCoverage.wig') as f:
        f.readline()
        f.readline()
        for num, line in enumerate(f):
            if int(line.rstrip()) < 20:
                seq[num] = 'n'
    seq = ''.join(seq)
    seq = seq.strip('n')
    with open(args.working_dir + '/%s.fasta' % args.sample, 'w') as o:
        o.write(">%s\n" % args.sample)
        for i in range(0, len(seq), 80):
            o.write(seq[i:i+80] + '\n')
    subprocess.Popen("prokka --force --cpus %s --outdir %s/prokka --prefix %s --kingdom Viruses --proteins %s/db/COVID.gbk  %s/%s.fasta "
                     % (args.threads, args.working_dir, args.sample, repo_dir, args.working_dir, args.sample), shell=True).wait()
    if bp / virus_length > 60:
        downsample = 60 / (bp / virus_length)
        subprocess.Popen("seqtk sample %s/reads.1.fq.gz %f | gzip > %s/reads.canu.1.fq.gz"
                         % (args.working_dir, downsample, args.working_dir), shell=True).wait()
        canu_reads_1 = "%s/reads.canu.1.fq.gz" % args.working_dir
    else:
        canu_reads_1 = "%s/reads.1.fq.gz" % args.working_dir
    subprocess.Popen("seqkit rmdup -s %s > %s/rmdup.fastq" % (canu_reads_1, args.working_dir), shell=True).wait()
    subprocess.Popen("canu -d %s/canu -pacbio-corrected %s/rmdup.fastq -p canu genomeSize=%d useGrid=false minOverlapLength=250"
                    % (args.working_dir, args.working_dir, virus_length), shell=True).wait()
    subprocess.Popen("minimap2 -t %s -ax map-pb %s/canu/canu.contigs.fasta %s | samtools view -b | samtools sort -@ %s -o %s/assembly.bam -"
                     " && samtools index %s/assembly.bam"
                     % (args.threads, args.working_dir, pilon_read_1, args.threads, args.working_dir, args.working_dir), shell=True).wait()
    subprocess.Popen("pilon --fix bases --threads %s --mindepth 20 --genome %s/canu/canu.contigs.fasta --unpaired %s/assembly.bam --tracks --output %s/assembly_pilon"
                     % (args.threads, args.working_dir, args.working_dir, args.working_dir), shell=True).wait()


__version__ = "0.1.1"
parser = argparse.ArgumentParser(prog='COVID pipeline', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='pipeline for the assembly, mapping, base calling of viruses\n' \
                                            'Version: %s\n' 
                                            'License: GPLv3\n'
                                            'USAGE: python -r1 <read1.fastq.gz> -r2 <read2.fastq.gz> -o sample1' % __version__)



parser.add_argument('-r1', '--read_1', action='store', help='Illumina pair 1')
parser.add_argument('-r2', '--read_2', action='store', help='Illumina pair 2')
parser.add_argument('-p', '--ccs_reads', action='store', help='Pacbio CCS reads')
parser.add_argument('-o', '--working_dir', action='store', default="temp", help='working directory')
parser.add_argument('-t', '--threads', action='store', default="12", help='number of threads to use')
parser.add_argument('-s', '--sample', action='store', help='sample name')
parser.add_argument('-c', '--coverage_pilon', default=200, type=int, action='store', help='downsample to this coverage for pilon')
parser.add_argument('-v', '--version', action='store_true', help="print version and exit")


args = parser.parse_args()


if args.version:
    sys.stdout.write('Version %s' % __version__)
    sys.exit()


if not os.path.exists(args.working_dir):
    os.makedirs(args.working_dir)
if not args.ccs_reads is None:
    run_ccs(args)
elif not args.read_1 is None and not args.read_2 is None:
    run_illumina(args)
else:
    sys.exit("Need read_1 and read_2 set or ccs_reads set")
