import os
import sys
import subprocess
import matplotlib
import argparse
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import statistics

def create_plots(sample_folder, amplified, threads, read1_suffix, read2_suffix, krakendb):
    repo_dir = sys.path[0]
    qc_dir = sample_folder + '/QC'
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)
    count = 0
    pp = PdfPages('%s/quality_control.pdf' % qc_dir)

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
    fig.suptitle("all reads\nconsensus %d bp long\n%d Ns" % (length, ns), fontsize=10)
    with open(qc_dir + '/report.txt', 'w') as o:
        o.write('== all reads ==\n')
        o.write("length: %d\n" % length)
        o.write("Ns: %d\n" % ns)
        o.write("max coverage: %d\n" % max(cov1))
        o.write("median coverage: %d\n" % statistics.median(cov1))
    for ax in axs.flat:
        ax.set(xlabel='position', ylabel='coverage')

    for ax in axs.flat:
        ax.label_outer()
    plt.savefig(pp, dpi=300)
    subprocess.Popen("samtools mpileup -f %s/db/COVID.fa %s/pipeline/ref.bam > %s/pileup" % (repo_dir, sample_folder, qc_dir), shell=True).wait()
    with open("%s/pileup" % qc_dir) as f:
        count_dict = {'a': [], 't': [], 'c': [], 'g': [], 'I': [], 'D': []}
        positions = []
        for line in f:
            ref, pos, refbase, cov, seq, qual = line.split()
            seq = seq.lower()
            counts = {'a':0, 't':0, 'c':0, 'g':0, 'I':0, 'D':0}
            depth = 0
            seq = list(seq)
            getdel = False
            getins = False
            while seq != []:
                x = seq.pop(0)
                mod = None
                if x == '.' or x == ',':
                    mod = refbase
                elif x == '+':
                    getins = True
                    digit = ''
                elif x == '-':
                    getdel = True
                    digit = ''
                elif x.isdigit() and (getdel or getins):
                    digit += x
                elif getdel:
                    if digit == '':
                        pass
                    else:
                        for j in range(int(digit) -1):
                            seq.pop(0)
                    mod = 'D'
                    getdel = False
                elif getins:
                    if digit == '':
                        pass
                    else:
                        for j in range(int(digit) -1):
                            seq.pop(0)
                    mod = 'I'
                    getins = False
                elif x in ['a', 't', 'c', 'g']:
                    mod = x
                if not mod is None:
                    counts[mod] += 1
                    depth += 1
            if depth >= 10 and counts[refbase] /depth < 0.9:
                for i in counts:
                    count_dict[i].append(counts[i])
                positions.append(pos + ' (' + refbase + ')')
    fig, ax = plt.subplots()
    for i in count_dict:
        ax.bar(positions, count_dict[i], 0.5, label=i)

    ax.set_ylabel('Counts')
    ax.set_title('Variants from reference')
    ax.set_xlabel('Position (ref. base)')

    plt.xticks(rotation=90)
    ax.legend()
    fig.subplots_adjust(bottom=0.2)
    plt.savefig(pp, dpi=300)
    # by read coverage bit
    for suffix in os.listdir(sample_folder):
        read1 = None
        read2 = None
        if suffix.endswith('.fastq'):
            continue
        for i in os.listdir(os.path.join(sample_folder, suffix)):
            if i.endswith(read1_suffix):
                read1 = os.path.join(sample_folder, suffix, i)
            if i.endswith(read2_suffix):
                read2 = os.path.join(sample_folder, suffix, i)
        if not read1 is None and not read2 is None:
            count += 1
    if count == 1:
        for i in os.listdir(repo_dir + '/db'):
            if i.startswith("SARS-CoV-2") and i.endswith(".csv"):
                with open(os.path.join(repo_dir, 'db', i)) as f:
                    labels = []
                    vals = []
                    xnums = []
                    f.readline()
                    num = 1
                    for line in f:
                        if line.split()[0].endswith("LEFT"):
                            label = line.split()[0][:-5]
                            l = int(line.split()[7])
                        elif line.split()[0].endswith("RIGHT"):
                            r = int(line.split()[7])
                            vals.append(statistics.median(cov1[l-1:r]))
                            xnums.append(num)
                            labels.append(label)
                            num += 1
                fig, ax = plt.subplots()
                ax.barh(xnums, vals)
                plt.yticks(xnums, labels)
                ax.set(ylabel="median depth", xlabel="primer set")
                ax.tick_params(axis='both', which='minor', labelsize=8)
                fig.suptitle("Median depth for primer set\n%s." % i[:-4], fontsize=10)

                plt.savefig(pp, dpi=300)
    else:
        for suffix in os.listdir(sample_folder):
            read1 = None
            read2 = None
            for i in os.listdir(os.path.join(sample_folder, suffix)):
                if i.endswith(read1_suffix):
                    read1 = os.path.join(sample_folder, suffix, i)
                if i.endswith(read2_suffix):
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
                fig.suptitle(suffix, fontsize=10)
                for ax in axs.flat:
                    ax.set(xlabel='position', ylabel='coverage')

                for ax in axs.flat:
                    ax.label_outer()
                plt.savefig(pp, dpi=300)
                for i in os.listdir(repo_dir + '/db'):
                    if i.startswith("SARS-CoV-2") and i.endswith(".csv"):
                        with open(os.path.join(repo_dir, 'db', i)) as f:
                            labels = []
                            vals = []
                            xnums = []
                            f.readline()
                            num = 1
                            for line in f:
                                if line.split()[0].endswith("LEFT"):
                                    label = line.split()[0][:-5]
                                    l = int(line.split()[7])
                                elif line.split()[0].endswith("RIGHT"):
                                    r = int(line.split()[7])
                                    vals.append(statistics.median(cov1[l-1:r]))
                                    xnums.append(num)
                                    labels.append(label)
                                    num += 1
                        fig, ax = plt.subplots()
                        ax.barh(xnums, vals)
                        plt.yticks(xnums, labels)
                        ax.tick_params(axis='both', which='minor', labelsize=8)
                        ax.set(ylabel="median depth", xlabel="primer set")
                        fig.suptitle("Median depth for primer set\n %s \nand reads %s." % (i[:-4], suffix), fontsize=10)
                        plt.savefig(pp, dpi=300)
                while True:
                    if cov1[0] >= 10:
                        break
                    cov1.pop(0)
                    if cov1 == []:
                        break
                gotzero = False
                while True:
                    count = 0
                    for i in cov1[-10:]:
                        if i < 10:
                            count += 1
                    if count < 3 and gotzero:
                        break
                    x = cov1.pop()
                    if x == 0:
                        gotzero = True
                    if cov1 == []:
                        break
                ns = 0
                for i in cov1:
                    if i <10:
                        ns += 1
                with open(qc_dir + '/report.txt', 'a') as o:
                    o.write("== %s ==\n" % suffix)
                    o.write("length: %d\n" % len(cov1))
                    o.write("Ns: %d\n" % ns)
                    o.write("max coverage: %d\n" % max(cov1))
                    o.write("median coverage: %d\n" % statistics.median(cov1))




    subprocess.Popen("samtools bam2fq -f 4 -@ %s %s/pipeline/ref.bam | gzip > %s/kraken_input.fastq.gz" % (threads, sample_folder, qc_dir), shell=True).wait()
    subprocess.Popen("kraken2 --db %s --quick --report %s/kraken_report.out --threads %s --output %s/kraken %s/kraken_input.fastq.gz" %
                     (krakendb, qc_dir, threads, qc_dir, qc_dir), shell=True).wait()
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
    fig.suptitle("kraken assignment", fontsize=10)
    plt.savefig(pp, dpi=300)
    pp.close()


__version__ = "0.1.1"
parser = argparse.ArgumentParser(prog='COVID pipeline QC', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='QC for the assembly, mapping, base calling of viruses\n' \
                                            'Version: %s\n' 
                                            'License: GPLv3\n'
                                            'USAGE: python run_QC.py -i sample1' % __version__)


parser.add_argument('-i', '--sample_folder', action='store', help='Sample folder created by process_run.py')
parser.add_argument('-a', '--not_amplified', action='store_true', help="Skip cutadapt")
parser.add_argument('-t', '--threads', action='store', default="12", help='number of threads to use')
parser.add_argument('-r1', '--read1_suffix', action='store', default="_1.fastq.gz", help='suffix for finding read 1')
parser.add_argument('-r2', '--read2_suffix', action='store', default="_2.fastq.gz", help='suffix for finding read 2')
parser.add_argument('-kdb', '--kraken_db', action='store', default="/sc/arion/projects/vanbah01b/COVID/db/minikraken2_v2_8GB_201904_UPDATE", help='location of kraken database')



args = parser.parse_args()

create_plots(args.sample_folder, not args.not_amplified, args.threads, args.read1_suffix, args.read2_suffix, args.kraken_db)
