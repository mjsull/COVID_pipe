import os
import sys
import subprocess
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import statistics

def create_plots(sample_folder, amplified, threads=12):
    repo_dir = sys.path[0]
    qc_dir = sample_folder + '/QC'
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)
    count = 0
    pp = PdfPages('multipage.pdf')

    # main coverage bit
    subprocess.Popen("samtools depth -aa %s/pipeline/ref.bam > %s/coverage.txt" % (sample_folder, qc_dir),
                     shell=True).wait()
    with open("%s/coverage.txt" % qc_dir) as f:
        cov1, cov2 = [], []
        for line in f:
            cov = int(line.split()[2])
            cov1.append(cov)
            cov2.append(min([100, cov]))
    fig, axs = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0})
    axs[0].plot(cov1)
    axs[1].plot(cov2)
    sample = os.path.basename(sample_folder.rstrip('/'))
    ns = 0
    length = 0
    with open(sample_folder + '/pipeline/' + sample + '.fasta') as f:
        for line in f:
            if not line.startswith('>'):
                length += len(line.rstrip())
                ns += line.lower().count('n')
    fig.suptitle("all reads\nconsenus %d bp\n%d Ns" % (length, ns), fontsize=16)
    for ax in axs.flat:
        ax.set(xlabel='position', ylabel='coverage')

    for ax in axs.flat:
        ax.label_outer()
    plt.savefig(pp, dpi=300)


    # by read coverage bit
    for suffix in os.listdir(sample_folder):
        read1 = None
        read2 = None
        if suffix.endswith('.fastq'):
            continue
        for i in os.listdir(os.path.join(sample_folder, suffix)):
            if i.endswith("R1_001.fastq.gz"):
                read1 = os.path.join(sample_folder, suffix, i)
            if i.endswith("R2_001.fastq.gz"):
                read2 = os.path.join(sample_folder, suffix, i)
        if not read1 is None and not read2 is None:
            count += 1
    if count == 1:
        for i in os.listdir(repo_dir + '/'):
            if i.startswith("SARS-CoV-2") and i.endswith(".tsv"):
                with open(os.path.join(repo_dir, i)) as f:
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
        fig, ax = plt.subplots()
        ax.bar(labels, vals, orientation="horizontal")
        fig.suptitle("Median coverage for primer set %s." % i[:-4], fontsize=16)
        plt.savefig(pp, dpi=300)
    else:
        for suffix in os.listdir(sample_folder):
            read1 = None
            read2 = None
            for i in os.listdir(os.path.join(sample_folder, suffix)):
                if i.endswith("R1_001.fastq.gz"):
                    read1 = os.path.join(sample_folder, suffix, i)
                if i.endswith("R2_001.fastq.gz"):
                    read2 = os.path.join(sample_folder, suffix, i)
            if not read1 is None and not read2 is None:
                if amplified:
                    subprocess.Popen(
                        "cutadapt -j %s -g file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -a file:%s/db/SARS-CoV-2_primers_3prime_NI.fa"
                        " -G file:%s/db/SARS-CoV-2_primers_5prime_NI.fa -A file:%s/db/SARS-CoV-2_primers_3prime_NI.fa "
                        "-o %s/reads.1.fq.gz -p %s/reads.2.fq.gz %s %s > %s/cutadapt.log"
                        % (threads, repo_dir, repo_dir, repo_dir, repo_dir, qc_dir, qc_dir,
                           read1, read2, qc_dir), shell=True).wait()
                    pilon_read_1 = "%s/reads.1.fq.gz" % qc_dir
                    pilon_read_2 = "%s/reads.2.fq.gz" % qc_dir
                else:
                    pilon_read_1 = read1
                    pilon_read_2 = read2
                subprocess.Popen(
                    "minimap2 -t %s -ax sr %s/db/COVID.fa %s %s | samtools view -b | samtools sort -@ %s -o %s/ref.bam -"
                    " && samtools index %s/ref.bam && samtools depth -aa %s/ref.bam > %s/coverage.txt"
                    % (threads, repo_dir, pilon_read_1, pilon_read_2, threads, qc_dir, qc_dir, qc_dir, qc_dir),
                    shell=True).wait()
                with open("%s/coverage.txt" % qc_dir) as f:
                    cov1, cov2 = [], []
                    for line in f:
                        cov = int(line.split()[2])
                        cov1.append(cov)
                        cov2.append(min([100, cov]))
                fig, axs = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0})
                axs[0].plot(cov1)
                axs[1].plot(cov2)
                fig.suptitle(suffix, fontsize=16)
                for ax in axs.flat:
                    ax.set(xlabel='position', ylabel='coverage')

                for ax in axs.flat:
                    ax.label_outer()
                plt.savefig(pp, dpi=300)
                for i in os.listdir(repo_dir + '/'):
                    if i.startswith("SARS-CoV-2") and i.endswith(".tsv"):
                        with open(os.path.join(repo_dir, i)) as f:
                            labels = []
                            vals = []
                            f.readline()
                            for line in f:
                                if line.split()[0].endswith("LEFT"):
                                    label = line.split()[0][:-5]
                                    l = int(line.split()[7])
                                elif line.split()[0].endswith("RIGHT"):
                                    r = int(line.split()[7])
                                    vals.append(statistics.median(cov1[l - 1:r]))
                                    labels.append(label)
                fig, ax = plt.subplots()
                ax.bar(labels, vals, orientation="horizontal")
                fig.suptitle("Median coverage for primer set %s with %s reads." % (i[:-4], suffix), fontsize=16)
                plt.savefig(pp, dpi=300)
    subprocess.Popen("samtools bam2fq -f 4 -@ %d %s/pipeline/ref.bam | gzip > %s/kraken_input.fastq.gz" % (threads, sample_folder, qc_dir), shell=True).wait()
    subprocess.Popen("kraken2 --db /sc/arion/projects/vanbah01b/COVID/db/minikraken2_v2_8GB_201904_UPDATE"
                     " --quick --report %s/kraken_report.out --threads %d --output %s/kraken %s/kraken_input.fastq.gz" % (qc_dir, threads, qc_dir, qc_dir), shell=True).wait()
    subprocess.Popen("samtools flagstat  %s/pipeline/ref.bam > %s/refbam.flagstat" % (sample_folder, qc_dir), shell=True).wait()
    with open("%s/refbam.flagstat" % qc_dir) as f:
        total_reads = int(f.readline().split()[0])
        f.readline()
        f.readline()
        f.readline()
        mapped_reads = int(f.readline().split()[0])
    chart = []
    chart.append((mapped_reads, "SARS-CoV-2"))
    other = 0
    with open("%s/kraken_report.out" % qc_dir) as f:
        for line in f:
            read_count = int(line.split()[2])
            classification = line.split()[5]
            if read_count / total_reads > 0.01:
                chart.append((read_count, classification))
            else:
                other += read_count
    chart.append((other, "other"))
    fig, ax = plt.subplots()
    sizes = []
    labels = []
    explode = []
    for i in chart:
        sizes.append(i[0]/total_reads*100)
        labels.append(i[1])
        explode.append(0)
    explode[0] = 0.1
    ax.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%', startangle=90)
    ax.axis('equal')
    fig.suptitle("kraken assignment", fontsize=16)
    plt.savefig(pp, dpi=300)






if len(sys.argv) == 3:
    create_plots(sys.argv[1], False)
else:
    create_plots(sys.argv[1], True)