import sys, os, argparse, subprocess
import pysam

def run_variant_analysis(args):
    sample_folder = args.sample_folder
    sample = os.path.basename(sample_folder.rstrip('/'))
    threads = args.threads
    outdir = sample_folder + '/variants'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    read_1 = "%s/pipeline/combined.1.fastq.gz" % sample_folder
    read_2 = "%s/pipeline/combined.2.fastq.gz" % sample_folder
    fasta = sample_folder + '/pipeline/' + sample + '.fasta'
    subprocess.Popen("minimap2 -t %s -ax sr %s %s %s | samtools view -b | samtools sort -@ %s -o %s/ref.bam -"
                     " && samtools index %s/ref.bam"
                     % (threads, fasta, read_1, read_2, args.threads, outdir, outdir), shell=True).wait()
    samfile = pysam.AlignmentFile("%s/ref.bam" % outdir, 'rb')
    outsam = pysam.AlignmentFile("%s/ref2.bam" % outdir, "wb", template=samfile)
    stored_reads_1 = {}
    stored_reads_2 = {}
    for read in samfile.fetch(until_eof=True):
        if not read.cigartuples is None:
            getit = True
            for i in read.cigartuples:
                if i[0] == 4 or i[0] == 5:
                    getit = False
                    break
            if not read.is_supplementary and not read.is_secondary and getit:
                if read.is_read1 and read.query_name in stored_reads_2:
                    outsam.write(stored_reads_2[read.query_name])
                    del stored_reads_2[read.query_name]
                    outsam.write(read)
                elif read.is_read1:
                    stored_reads_1[read.query_name] = read
                elif read.is_read2 and read.query_name in stored_reads_1:
                    outsam.write(stored_reads_1[read.query_name])
                    del stored_reads_1[read.query_name]
                    outsam.write(read)
                elif read.is_read2:
                    stored_reads_2[read.query_name] = read
    outsam.close()
    subprocess.Popen("samtools sort -o %s/ref3.bam %s/ref2.bam && samtools index %s/ref3.bam && samtools mpileup -f %s %s/ref3.bam > %s/pileup" % (outdir, outdir, outdir, fasta, outdir, outdir), shell=True).wait()
    modlist = ['a', 't', 'c', 'g', 'n', 'I', 'D']
    outseq = ''
    with open("%s/pileup" % outdir) as f, open("%s/variable_bases.tsv" % outdir, "w") as out:
        out.write("reference\tposition\treference_base\tcoverage\treference_base_fraction\tflagged\ta\tt\tc\tg\tn\tinsertion\tdeletion\n")
        for line in f:
            ref, pos, refbase, cov, seq, qual = line.split()
            refbase = refbase.lower()
            seq = seq.lower()
            counts = {'a':0, 't':0, 'c':0, 'g':0, 'n':0, 'I':0, 'D':0}
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
                    for j in range(int(digit) -1):
                        seq.pop(0)
                    mod = 'D'
                    getdel = False
                elif getins:
                    for j in range(int(digit) -1):
                        seq.pop(0)
                    mod = 'I'
                    getins = False
                elif x == '^':
                    seq.pop(0)
                elif x in ['a', 't', 'c', 'g']:
                    mod = x
                elif x in ['$', '*']:
                    pass
                else:
                    sys.exit("Base not recognised " + x)
                if not mod is None:
                    counts[mod] += 1
                    depth += 1
            try:
                fraction = counts[refbase]/depth
            except ZeroDivisionError:
                fraction = None
            if fraction is None:
                flagged = "no_cov"
                outseq += refbase
            elif fraction < args.min_ratio:
                flagged = "FLAGGED"
                outseq += 'n'
            else:
                flagged = "OK"
                outseq += refbase
            outlist = [ref, pos, refbase, cov, fraction, flagged]
            for i in modlist:
                outlist.append(counts[i])
            out.write('\t'.join(map(str, outlist)) + '\n')
    with open("%s/%s.final.fna" % (outdir, sample), "w") as out:
        out.write(">%s\n" % sample)
        for i in range(0, len(outseq), 80):
            out.write(outseq[i:i+80] + '\n')




__version__ = "0.1.1"
parser = argparse.ArgumentParser(prog='COVID variant pipeline', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='QC for the assembly, mapping, base calling of viruses\n' \
                                            'Version: %s\n' 
                                            'License: GPLv3\n'
                                            'USAGE: python run_QC.py -i sample1' % __version__)


parser.add_argument('-i', '--sample_folder', action='store', help='Sample folder created by process_run.py')
parser.add_argument('-t', '--threads', action='store', default="12", help='number of threads to use')
parser.add_argument('-m', '--min_ratio', action='store', default=0.85, help='number of threads to use')



args = parser.parse_args()
run_variant_analysis(args)