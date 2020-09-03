import argparse
import sys
import os
import subprocess


def get_variants(args):
    seq_dict = {}
    ref_name = args.reference_name
    out_dict = {}
    with open(args.multiple_alignment) as f:
        for line in f:
            if line.startswith(">"):
                name = line.split()[0][1:]
                seq_dict[name] = ''
                out_dict[name] = {}
            else:
                seq_dict[name] += line.rstrip().lower()
    with open(args.lineage_file) as f:
        lineage_dict = {}
        lineage_count = {}
        for line in f:
            name, lineage = line.split(',')[:2]
            lineage_dict[name] = lineage
            if not lineage in lineage_count:
                lineage_count[lineage] = 0
            lineage_count[lineage] += 1
    count = 1
    for i in range(0, len(seq_dict[ref_name])):
        if seq_dict[ref_name][i] == '-':
            pass
        else:
            seqlist = []
            for j in seq_dict:
                seqlist.append(seq_dict[j][i])
            bases_at_position = []
            if seqlist.count('a') >= args.min_base:
                bases_at_position.append('a')
            if seqlist.count('t') >= args.min_base:
                bases_at_position.append('t')
            if seqlist.count('c') >= args.min_base:
                bases_at_position.append('c')
            if seqlist.count('g') >= args.min_base:
                bases_at_position.append('g')
            if len(bases_at_position) >= 2:
                out_dict[count] = {}
                for base in bases_at_position:
                    count_dict = {}
                    for j in seq_dict:
                        if seq_dict[j][i] == base:
                            lineage = lineage_dict[j]
                            if not lineage in count_dict:
                                count_dict[lineage] = 0
                            count_dict[lineage] += 1
                    outstring = []
                    for j in count_dict:
                        outstring.append(j + ':%.2f' % (count_dict[j]/lineage_count[j]*100))
                    out_dict[count][base] = ';'.join(outstring)
            count += 1
    return out_dict, seq_dict[args.reference_name]

def get_profile_pileup(args, var_dict):
    major_allele = {}
    minor_allele = {}
    with open(args.pileup) as f, open(args.out_file, 'w') as o:
        for line in f:
            ref, pos, refbase, cov, seq, qual = line.split()
            refbase = refbase.lower()
            seq = seq.lower()
            counts = {'a':0, 't':0, 'c':0, 'g':0, 'n':0, 'I':0, 'D':0}
            depth = 0
            seq = list(seq)
            pos = int(pos)
            minor_allele[pos] = set()
            major_allele[pos] = set()
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
                elif x == '-' and not getdel:
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
            if depth >= args.min_depth:
                if pos in var_dict:
                    varout = [str(pos), str(depth)]
                    bases = ['a', 't', 'c', 'g']
                    bases.sort(key=lambda x: counts[x], reverse=True)
                    for i in bases:
                        if counts[i] /depth > args.min_fraction and counts[i] > args.min_variant_count and i in var_dict[pos]:
                            varout += [i, var_dict[pos][i], str(counts[i])]
                    if len(varout) > 5:
                        o.write("\t".join(varout) + "\n")

    return (major_allele, minor_allele)

__version__ = "0.0.1"

parser = argparse.ArgumentParser(prog='contamination_analysis.py', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='Show where hetrozygous bases show up on a tree.\n' \
                                            'Version: %s\n' 
                                            'License: GPLv3\n'
                                            'USAGE: python -p output.pileup -m multiple_alignment -r reference_name -t tree.nwk' % __version__)

parser.add_argument('-i', '--input_dir', action='store', help='automatically run on directory')
parser.add_argument('-o', '--out_file', action='store', help='automatically run on directory')
parser.add_argument('-m', '--multiple_alignment', action='store', help='number of threads to use')
parser.add_argument('-l', '--lineage_file', action='store', help='lineage file')
parser.add_argument('-r', '--reference_name', action='store', help='name of reference in multiple alignment')
parser.add_argument('-p', '--pileup', action='store', help='pileup of reads')
parser.add_argument('-f', '--min_fraction', action='store', default=0.05, help='minimum fraction of alternate to flag')
parser.add_argument('-c', '--min_variant_count', action='store', default=2, help='minimum depth to consider site')
parser.add_argument('-d', '--min_depth', action='store', default=50, help='minimum depth to consider site')
parser.add_argument('-b', '--min_base', action='store', default=2, help='minimum isolates base included in to include variant')



args = parser.parse_args()

if not args.input_dir is None:
    repo_dir = sys.path[0]
    working_dir = args.input_dir + "/contamination_analysis"
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    args.out_file = working_dir + "/report.tsv"
    args.pileup = working_dir + "/pileup"
    read1 = os.path.join(args.input_dir, "pipeline", "combined.1.fastq.gz")
    read2 = os.path.join(args.input_dir, "pipeline", "combined.2.fastq.gz")
    ref_fasta = working_dir + '/ref.fa'
    args.multiple_alignment = os.path.join(repo_dir, 'db', 'multiple_alignment.fa')
    args.lineage_file = os.path.join(repo_dir, 'db', 'lineage_report.csv')
    args.reference_name = "Sequence_Identifier:VS_7528.aid_7029.SARS-CoV-2|Sample_Identifier:VS_7528|Sample_Name:PV13343|Segment:SARS-CoV-2"
    var_dict, ref_seq = get_variants(args)
    with open(ref_fasta, 'w') as o:
        o.write(">" + args.reference_name + "\n" + ref_seq.replace('-', '') + "\n")
    subprocess.Popen("minimap2 -ax sr -t 8 %s %s %s | samtools view -b | samtools sort -o %s/var_tree.bam -" %
                     (ref_fasta, read1, read2, working_dir), shell=True).wait()
    subprocess.Popen("samtools mpileup -f %s %s/var_tree.bam > %s" % (ref_fasta, working_dir, args.pileup), shell=True).wait()
    get_profile_pileup(args, var_dict)

