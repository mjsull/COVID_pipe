#!/usr/bin/env python3

import argparse
import subprocess
import sys, os


repo_dir = sys.path[0]
virus_length = 29760

def run_illumina(args):
    read1s = []
    read2s = []
    working_dir = os.path.join(args.sample_folder, "pipeline")
    sample = os.path.basename(args.sample_folder.rstrip('/'))
    ion_reads = None
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    for suffix in os.listdir(args.sample_folder):
        if suffix.endswith('Z.fastq'):
            ion_reads = os.path.join(args.sample_folder, suffix)
        else:
            for i in os.listdir(os.path.join(args.sample_folder, suffix)):
                if i.endswith("R1_001.fastq.gz"):
                    read1s.append(os.path.join(args.sample_folder, suffix, i))
                if i.endswith("R2_001.fastq.gz"):
                    read2s.append(os.path.join(args.sample_folder, suffix, i))
    if not ion_reads is None:
        pass
    elif len(read1s) == 0 or len(read2s) == 0:
        sys.exit("Couldn't find any reads in sample folder")
    else:
        subprocess.Popen("cat %s  > %s/combined.1.fastq.gz" % (' '.join(read1s), working_dir), shell=True).wait()
        subprocess.Popen("cat %s  > %s/combined.2.fastq.gz" % (' '.join(read2s), working_dir), shell=True).wait()
    if not ion_reads is None:
        print("cutadapt -j %s -g file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -a file:%s/db/SARS-CoV-2_primers_3prime_NI.fa -o %s/reads.1.fq.gz %s > %s/cutadapt.1.log" % (args.threads, repo_dir, repo_dir, working_dir, ion_reads, working_dir))
        subprocess.Popen("cutadapt -j %s -g file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -a file:%s/db/SARS-CoV-2_primers_3prime_NI.fa -o %s/reads.1.fq.gz %s > %s/cutadapt.1.log" % (args.threads, repo_dir, repo_dir, working_dir, ion_reads, working_dir), shell=True).wait()
        pilon_read_1 = "%s/reads.1.fq.gz" % working_dir
        pilon_read_2 = ""
    elif not args.not_amplified:
        subprocess.Popen("cutadapt -j %s -g file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -a file:%s/db/SARS-CoV-2_primers_3prime_NI.fa"
                     " -G file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -A file:%s/db/SARS-CoV-2_primers_3prime_NI.fa "
                     "-o %s/reads.1.fq.gz -p %s/reads.2.fq.gz %s/combined.1.fastq.gz %s/combined.2.fastq.gz > %s/cutadapt.1.log"
                     % (args.threads, repo_dir, repo_dir, repo_dir, repo_dir, working_dir, working_dir,
                        working_dir, working_dir, working_dir), shell=True).wait()
        pilon_read_1 = "%s/reads.1.fq.gz" % working_dir
        pilon_read_2 = "%s/reads.2.fq.gz" % working_dir
    else:
        pilon_read_1 = "%s/combined.1.fastq.gz" % working_dir
        pilon_read_2 = "%s/combined.2.fastq.gz" % working_dir
    subprocess.Popen("minimap2 -t %s -ax sr %s/db/COVID.fa %s %s | samtools view -b | samtools sort -@ %s -o %s/ref.bam -"
                     " && samtools index %s/ref.bam"
                     % (args.threads, repo_dir, pilon_read_1, pilon_read_2, args.threads, working_dir, working_dir), shell=True).wait()
    subprocess.Popen("pilon --fix bases --changes --vcf --threads %s --mindepth 10 --genome %s/db/COVID.fa --frags %s/ref.bam --tracks --output %s/pilon"
                     % (args.threads, repo_dir, working_dir, working_dir), shell=True).wait()
    with open(working_dir + '/pilon.fasta') as f:
        seq = ''
        for line in f:
            if not line.startswith('>'):
                seq += line.rstrip()
    seq = list(seq)
    with open(working_dir + '/pilonCoverage.wig') as f:
        f.readline()
        f.readline()
        for num, line in enumerate(f):
            if int(line.rstrip()) < 10:
                seq[num] = 'n'
    seq = ''.join(seq)
    if seq.endswith('a'):
        seq = seq.rstrip('a')
    seq = seq.strip('n')
    while True:
        if len(seq) < 10:
            break
        if seq[-10:].count('n') < 3:
            break
        seq = seq[:-1]
        seq = seq.rstrip('n')
    with open(working_dir + '/%s.fasta' % sample, 'w') as o:
        o.write(">%s\n" % sample)
        for i in range(0, len(seq), 80):
            o.write(seq[i:i + 80] + '\n')
    subprocess.Popen(
        "prokka --force --cpus %s --outdir %s/prokka --prefix %s --kingdom Viruses --proteins %s/db/COVID.gbk  %s/%s.fasta "
        % (args.threads, working_dir, sample, repo_dir, working_dir, sample), shell=True).wait()
    if not args.not_amplified and ion_reads is None:
        subprocess.Popen("shovill --outdir %s/shovill --R1 %s --R2 %s --gsize %d --cpus %s"
                     % (working_dir, pilon_read_1, pilon_read_2, virus_length, args.threads), shell=True).wait()



def run_ccs(args):
    subprocess.Popen("cutadapt -j %s -g file:%s/db/SARS-CoV-2_primers_5prime_anchored.fa -a "
                     "file:%s/db/SARS-CoV-2_primers_3prime_anchored.fa -o %s/reads.1.fq.gz %s > %s/cutadapt.1.log"
                     % (args.threads, repo_dir, repo_dir, working_dir, args.ccs_reads, working_dir), shell=True).wait()
    with open(working_dir + '/cutadapt.1.log') as f:
        for line in f:
            if line.startswith("Total written (filtered):"):
                bp = int(line.split()[3].replace(',', ''))
                break
    if bp / virus_length > args.coverage_pilon:
        downsample = args.coverage_pilon / (bp / virus_length)
        subprocess.Popen("seqtk sample %s/reads.1.fq.gz %f | gzip > %s/reads.pilon.1.fq.gz"
                         % (working_dir, downsample, working_dir), shell=True).wait()
        pilon_read_1 = "%s/reads.pilon.1.fq.gz" % working_dir
    else:
        pilon_read_1 = "%s/reads.1.fq.gz" % working_dir
    subprocess.Popen("minimap2 -t %s -ax map-pb %s/db/COVID.fa %s | samtools view -b | samtools sort -@ %s -o %s/ref.bam -"
                     " && samtools index %s/ref.bam"
                     % (args.threads, repo_dir, pilon_read_1, args.threads, working_dir, working_dir), shell=True).wait()
    subprocess.Popen("pilon --fix bases --threads %s --mindepth 20 --genome %s/db/COVID.fa --unpaired %s/ref.bam --tracks --output %s/pilon"
                     % (args.threads, repo_dir, working_dir, working_dir), shell=True).wait()
    with open(working_dir + '/pilon.fasta') as f:
        seq = ''
        for line in f:
            if not line.startswith('>'):
                seq += line.rstrip()
    seq = list(seq)
    with open(working_dir + '/pilonCoverage.wig') as f:
        f.readline()
        f.readline()
        for num, line in enumerate(f):
            if int(line.rstrip()) < 20:
                seq[num] = 'n'
    seq = ''.join(seq)
    seq = seq.strip('n')
    with open(working_dir + '/%s.fasta' % args.sample, 'w') as o:
        o.write(">%s\n" % args.sample)
        for i in range(0, len(seq), 80):
            o.write(seq[i:i+80] + '\n')
    subprocess.Popen("prokka --force --cpus %s --outdir %s/prokka --prefix %s --kingdom Viruses --proteins %s/db/COVID.gbk  %s/%s.fasta "
                     % (args.threads, working_dir, args.sample, repo_dir, working_dir, args.sample), shell=True).wait()
    if bp / virus_length > 60:
        downsample = 60 / (bp / virus_length)
        subprocess.Popen("seqtk sample %s/reads.1.fq.gz %f | gzip > %s/reads.canu.1.fq.gz"
                         % (working_dir, downsample, working_dir), shell=True).wait()
        canu_reads_1 = "%s/reads.canu.1.fq.gz" % working_dir
    else:
        canu_reads_1 = "%s/reads.1.fq.gz" % working_dir
    subprocess.Popen("seqkit rmdup -s %s > %s/rmdup.fastq" % (canu_reads_1, working_dir), shell=True).wait()
    subprocess.Popen("canu -d %s/canu -pacbio-corrected %s/rmdup.fastq -p canu genomeSize=%d useGrid=false minOverlapLength=250"
                    % (working_dir, working_dir, virus_length), shell=True).wait()
    subprocess.Popen("minimap2 -t %s -ax map-pb %s/canu/canu.contigs.fasta %s | samtools view -b | samtools sort -@ %s -o %s/assembly.bam -"
                     " && samtools index %s/assembly.bam"
                     % (args.threads, working_dir, pilon_read_1, args.threads, working_dir, working_dir), shell=True).wait()
    subprocess.Popen("pilon --fix bases --threads %s --mindepth 20 --genome %s/canu/canu.contigs.fasta --unpaired %s/assembly.bam --tracks --output %s/assembly_pilon"
                     % (args.threads, working_dir, working_dir, working_dir), shell=True).wait()


__version__ = "0.1.1"
parser = argparse.ArgumentParser(prog='COVID pipeline', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='pipeline for the assembly, mapping, base calling of viruses\n' \
                                            'Version: %s\n' 
                                            'License: GPLv3\n'
                                            'USAGE: python -r1 <read1.fastq.gz> -r2 <read2.fastq.gz> -o sample1' % __version__)



parser.add_argument('-i', '--sample_folder', action='store', help='Sample folder created by process_run.py')
parser.add_argument('-p', '--ccs_reads', action='store', help='Pacbio CCS reads')
parser.add_argument('-o', '--working_dir', action='store', default="temp", help='working directory')
parser.add_argument('-t', '--threads', action='store', default="12", help='number of threads to use')
parser.add_argument('-s', '--sample', action='store', help='sample name')
parser.add_argument('-c', '--coverage_pilon', default=200, type=int, action='store', help='downsample to this coverage for pilon')
parser.add_argument('-v', '--version', action='store_true', help="print version and exit")
parser.add_argument('-a', '--not_amplified', action='store_true', help="Skip cutadapt and assembly")


args = parser.parse_args()


if args.version:
    sys.stdout.write('Version %s' % __version__)
    sys.exit()


if not args.ccs_reads is None:
    run_ccs(args)
elif not args.sample_folder is None:
    run_illumina(args)
else:
    sys.exit("Need read_1 and read_2 set or ccs_reads set")