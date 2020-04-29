import sys
import subprocess
import os
import statistics

read1_suffix = "_1.fastq.gz"
read2_suffix = "_2.fastq.gz"


qc_dir = sys.argv[2]
try:
    os.makedirs(qc_dir)
except:
    pass

repo_dir = sys.path[0]
threads = 24
project_folder = sys.argv[1]


first = True
with open(sys.argv[3] + '.15.tsv', 'w') as o1, open(sys.argv[3] + '.2.tsv', 'w') as o2:
    for qq in os.listdir(project_folder):
        sample_folder = os.path.join(project_folder, qq)
        if sample_folder == "H2z0" or not os.path.isdir(sample_folder):
            continue
        read1_15 = None
        read2_15 = None
        read1_2 = None
        read2_2 = None
        for suffix in os.listdir(sample_folder):
            read1 = None
            read2 = None
            if not os.path.isdir(os.path.join(sample_folder, suffix)) or suffix == 'pipeline' or suffix == 'QC':
                continue
            for i in os.listdir(os.path.join(sample_folder, suffix)):
                if i.endswith(read1_suffix):
                    read1 = os.path.join(sample_folder, suffix, i)
                if i.endswith(read2_suffix):
                    read2 = os.path.join(sample_folder, suffix, i)
            if not read1 is None and not read2 is None and "5kb" in suffix.lower():
                read1_15 = read1
                read2_15 = read2
            elif not read1 is None and not read2 is None and "0kb" in suffix.lower():
                read1_2 = read1
                read2_2 = read2
        if not None in [read1_2, read2_2, read1_15, read2_15]:
            subprocess.Popen(
                    "cutadapt -j %s -g file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -a file:%s/db/SARS-CoV-2_primers_3prime_NI.fa"
                    " -G file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -A file:%s/db/SARS-CoV-2_primers_3prime_NI.fa "
                    "-o %s/reads.1.fq.gz -p %s/reads.2.fq.gz %s %s > %s/cutadapt.log"
                    % (threads, repo_dir, repo_dir, repo_dir, repo_dir, qc_dir, qc_dir,
                       read1_15, read2_15, qc_dir), shell=True).wait()
            pilon_read_1 = "%s/reads.1.fq.gz" % qc_dir
            pilon_read_2 = "%s/reads.2.fq.gz" % qc_dir
            subprocess.Popen(
                "minimap2 -t %s -ax sr %s/db/COVID.fa %s %s | samtools view -b | samtools sort -@ %s -o %s/ref.bam -"
                " && samtools index %s/ref.bam && samtools depth -aa %s/ref.bam > %s/coverage.txt"
                % (threads, repo_dir, pilon_read_1, pilon_read_2, threads, qc_dir, qc_dir, qc_dir, qc_dir),
                shell=True).wait()
            with open("%s/coverage.txt" % qc_dir) as f:
                cov1 = []
                for line in f:
                    cov = int(line.split()[2])
                    cov1.append(cov)
            with open(os.path.join(repo_dir, 'db/SARS-CoV-2_primers_1.5kb_set.csv')) as f:
                labels = []
                vals = []
                f.readline()
                for line in f:
                    if line.split()[0].endswith("LEFT"):
                        label = line.split()[0][:-5]
                        l = int(line.split()[7])
                    elif line.split()[0].endswith("RIGHT"):
                        r = int(line.split()[7])
                        vals.append(statistics.median(cov1[l-1:r]))
                        labels.append(label)
            if first:
                o1.write('sample\t' + '\t'.join(labels) + '\n')
            o1.write(qq + '\t' + '\t'.join(map(str, vals)) + '\n')
            subprocess.Popen(
                    "cutadapt -j %s -g file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -a file:%s/db/SARS-CoV-2_primers_3prime_NI.fa"
                    " -G file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -A file:%s/db/SARS-CoV-2_primers_3prime_NI.fa "
                    "-o %s/reads.1.fq.gz -p %s/reads.2.fq.gz %s %s > %s/cutadapt.log"
                    % (threads, repo_dir, repo_dir, repo_dir, repo_dir, qc_dir, qc_dir,
                       read1_2, read2_2, qc_dir), shell=True).wait()
            pilon_read_1 = "%s/reads.1.fq.gz" % qc_dir
            pilon_read_2 = "%s/reads.2.fq.gz" % qc_dir
            subprocess.Popen(
                "minimap2 -t %s -ax sr %s/db/COVID.fa %s %s | samtools view -b | samtools sort -@ %s -o %s/ref.bam -"
                " && samtools index %s/ref.bam && samtools depth -aa %s/ref.bam > %s/coverage.txt"
                % (threads, repo_dir, pilon_read_1, pilon_read_2, threads, qc_dir, qc_dir, qc_dir, qc_dir),
                shell=True).wait()
            with open("%s/coverage.txt" % qc_dir) as f:
                cov1 = []
                for line in f:
                    cov = int(line.split()[2])
                    cov1.append(cov)
            with open(os.path.join(repo_dir, 'db/SARS-CoV-2_primers_2kb_set.csv')) as f:
                labels = []
                vals = []
                f.readline()
                for line in f:
                    if line.split()[0].endswith("LEFT"):
                        label = line.split()[0][:-5]
                        l = int(line.split()[7])
                    elif line.split()[0].endswith("RIGHT"):
                        r = int(line.split()[7])
                        vals.append(statistics.median(cov1[l-1:r]))
                        labels.append(label)
            if first:
                first = False
                o2.write('sample\t' + '\t'.join(labels) + '\n')
            o2.write(qq + '\t' + '\t'.join(map(str, vals)) + '\n')




