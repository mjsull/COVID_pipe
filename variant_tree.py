import argparse, sys
from ete3 import Tree, RectFace, AttrFace, TextFace
from ete3 import NodeStyle
from ete3 import TreeStyle
import os
import subprocess

def colorstr(rgb): return "#%02x%02x%02x" % (rgb[0],rgb[1],rgb[2])


def hsl_to_rgb(h, s, l):
    c = (1 - abs(2*l - 1)) * s
    x = c * (1 - abs(h *1.0 / 60 % 2 - 1))
    m = l - c/2
    if h < 60:
        r, g, b = c + m, x + m, 0 + m
    elif h < 120:
        r, g, b = x + m, c+ m, 0 + m
    elif h < 180:
        r, g, b = 0 + m, c + m, x + m
    elif h < 240:
        r, g, b, = 0 + m, x + m, c + m
    elif h < 300:
        r, g, b, = x + m, 0 + m, c + m
    else:
        r, g, b, = c + m, 0 + m, x + m
    r = int(r * 255)
    g = int(g * 255)
    b = int(b * 255)
    return (r,g,b)



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
    count = 1
    for i in range(0, len(seq_dict[ref_name])):
        if seq_dict[ref_name] == '-':
            pass
        else:
            seqlist = []
            for j in seq_dict:
                seqlist.append(seq_dict[j][i])
            bases_at_position = 0
            if seqlist.count('a') >= args.min_base:
                bases_at_position += 1
            if seqlist.count('t') >= args.min_base:
                bases_at_position += 1
            if seqlist.count('c') >= args.min_base:
                bases_at_position += 1
            if seqlist.count('g') >= args.min_base:
                bases_at_position += 1
            if bases_at_position >= 2:
                for j in seq_dict:
                    if seq_dict[j][i] == 'a':
                        out_dict[j][count] = 'a'
                    elif seq_dict[j][i] == 't':
                        out_dict[j][count] = 't'
                    elif seq_dict[j][i] == 'c':
                        out_dict[j][count] = 'c'
                    elif seq_dict[j][i] == 'g':
                        out_dict[j][count] = 'g'
                    else:
                        out_dict[j][count] = 'x'
            count += 1
    return out_dict


def get_profile_pileup(args):
    major_allele = {}
    minor_allele = {}
    with open(args.pileup) as f:
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
            if depth >= args.min_depth:
                for i in ['a', 't', 'c', 'g']:
                    if counts[i] /depth > args.min_fraction:
                        if counts[i]/depth > 0.5:
                            major_allele[pos].add(i)
                        else:
                            minor_allele[pos].add(i)
    return (major_allele, minor_allele)





def draw_tree(variant_dict, profile_dict, the_tree, ref_name):
    major_alleles, minor_alleles = profile_dict
    t = Tree(the_tree, quoted_node_names=True)
    font_size = 8
    font_type = 'Heveltica'
    font_gap = 3
    font_buffer = 10
    ts = TreeStyle()
    position_list = None
    allele_count = {}
    max_major = 0
    max_minor = 0
    max_combo = 0
    max_major_minor = 0
    for n in t.iter_leaves():
        the_name = n.name
        if the_name == ref_name:
            t.set_outgroup(n)
        if position_list is None:
            position_list = list(variant_dict[the_name])
            position_list.sort()
            for num, i in enumerate(position_list):
                nameF = TextFace(font_gap * ' ' + str(i) + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True)
                nameF.rotation = -90
                ts.aligned_header.add_face(nameF, column=num)
        minor_allele, major_allele, missing = 0, 0, 0
        for num, i in enumerate(position_list):
            if i in major_alleles and variant_dict[the_name][i] in major_alleles[i]:
                s, l = 0.0, 0.1
                major_allele += 1
            elif i in minor_alleles and variant_dict[the_name][i] in minor_alleles[i]:
                s, l = 0.0, 0.5
                minor_allele += 1
            else:
                s, l = 0.0, 0.9
                missing += 1
            if variant_dict[the_name][i] == 'a':
                colour = colorstr(hsl_to_rgb(0, s, l))
                #colour = "#65dad0"
            elif variant_dict[the_name][i] == 't':
                colour = colorstr(hsl_to_rgb(140, s, l))
                #colour = "#daa3dc"
            elif variant_dict[the_name][i] == 'c':
                colour = colorstr(hsl_to_rgb(220, s, l))
                #colour = "#9bd686"
            elif variant_dict[the_name][i] == 'g':
                colour = colorstr(hsl_to_rgb(300, s, l))
                #colour = "#e1b86f"
            else:
                colour = "#ffffff"
            #colour = colorstr(hsl_to_rgb(h, s, l))
            n.add_face(RectFace(20, 20, colour, colour), column=num, position="aligned")
        allele_count[the_name] = (major_allele, minor_allele, missing)
        if major_allele > max_major:
            max_major = major_allele
        if minor_allele + major_allele >= max_combo:
            max_combo = minor_allele + major_allele
        if minor_allele >= max_minor:
            max_minor = minor_allele
            if major_allele > max_major_minor:
                max_major_minor = major_allele
    nst1 = NodeStyle()
    nst1["bgcolor"] = "LightSteelBlue"
    nst2 = NodeStyle()
    nst2["bgcolor"] = "Moccasin"
    with open(args.out_text, 'w') as o:
        for n in t.iter_leaves():
            the_name = n.name
            major_allele, minor_allele, missing = allele_count[the_name]
            if major_allele == max_major:
                o.write("max_major\t%s\t%d\t%d\t%d\n" % (the_name,major_allele,minor_allele, missing))
                n.img_style["bgcolor"] = "#cb6a49"
            elif minor_allele == max_minor and major_allele == max_major_minor:
                o.write("max_minor\t%s\t%d\t%d\t%d\n" % (the_name,major_allele,minor_allele, missing))
                n.img_style["bgcolor"] = "#7aa457"
            elif minor_allele + major_allele == max_combo:
                o.write("max_combo\t%s\t%d\t%d\t%d\n" % (the_name,major_allele,minor_allele, missing))
                n.img_style["bgcolor"] = "#a46cb7"
            n.add_face(TextFace('%d/%d/%d' % (major_allele, minor_allele, missing), fsize=font_size, ftype=font_type,
                     tight_text=True), column=len(position_list), position="aligned")



    ts.legend.add_face(
        TextFace('A', fsize=font_size, ftype=font_type,
                 tight_text=True), column=0)
    s = 0.5
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(0, s, 0.3)), colorstr(hsl_to_rgb(0, s, 0.3))), column=1)
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(0, s, 0.5)), colorstr(hsl_to_rgb(0, s, 0.5))), column=2)
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(0, s, 0.8)), colorstr(hsl_to_rgb(0, s, 0.8))), column=3)
    ts.legend.add_face(
        TextFace('T', fsize=font_size, ftype=font_type,
                 tight_text=True), column=0)
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(140, s, 0.3)), colorstr(hsl_to_rgb(140, s, 0.3))), column=1)
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(140, s, 0.5)), colorstr(hsl_to_rgb(140, s, 0.5))), column=2)
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(140, s, 0.8)), colorstr(hsl_to_rgb(140, s, 0.8))), column=3)
    ts.legend.add_face(
        TextFace('C', fsize=font_size, ftype=font_type,
                 tight_text=True), column=0)
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(220, s, 0.3)), colorstr(hsl_to_rgb(220, s, 0.3))), column=1)
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(220, s, 0.5)), colorstr(hsl_to_rgb(220, s, 0.5))), column=2)
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(220, s, 0.8)), colorstr(hsl_to_rgb(220, s, 0.8))), column=3)
    ts.legend.add_face(
        TextFace('G', fsize=font_size, ftype=font_type,
                 tight_text=True), column=0)
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(300, s, 0.3)), colorstr(hsl_to_rgb(300, s, 0.3))), column=1)
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(300, s, 0.5)), colorstr(hsl_to_rgb(300, s, 0.5))), column=2)
    ts.legend.add_face(RectFace(20, 20, colorstr(hsl_to_rgb(300, s, 0.8)), colorstr(hsl_to_rgb(300, s, 0.8))), column=3)
    ts.legend.add_face(
        TextFace('-', fsize=font_size, ftype=font_type,
                 tight_text=True), column=0)
    ts.legend.add_face(RectFace(20, 20, "#cccccc", "#cccccc"), column=1)
    t.render(args.out_file, w=210, units='mm', tree_style=ts)


__version__ = "0.0.1"





parser = argparse.ArgumentParser(prog='Variant_tree.py', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='Show where hetrozygous bases show up on a tree.\n' \
                                            'Version: %s\n' 
                                            'License: GPLv3\n'
                                            'USAGE: python -p output.pileup -m multiple_alignment -r reference_name -t tree.nwk' % __version__)

parser.add_argument('-i', '--input_dir', action='store', help='automatically run on directory')
parser.add_argument('-o', '--out_file', action='store', help='automatically run on directory')
parser.add_argument('-m', '--multiple_alignment', action='store', help='number of threads to use')
parser.add_argument('-r', '--reference_name', action='store', help='name of reference in multiple alignment')
parser.add_argument('-t', '--tree', action='store', help='number of threads to use')
parser.add_argument('-p', '--pileup', action='store', help='pileup of reads')
parser.add_argument('-f', '--min_fraction', action='store', default=0.05, help='minimum fraction of alternate to flag')
parser.add_argument('-d', '--min_depth', action='store', default=20, help='minimum depth to consider site')
parser.add_argument('-b', '--min_base', action='store', default=2, help='minimum isolates base included in to include variant')



args = parser.parse_args()

if not args.input_dir is None:
    repo_dir = sys.path[0]
    working_dir = args.input_dir + "/var_tree"
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    args.out_file = working_dir + "/tree.pdf"
    args.pileup = working_dir + "/pileup"
    args.multiple_alignment = repo_dir + '/db/ncov_MSA_2020-05-19_237.fasta'
    args.out_text = working_dir + "/description.txt"
    ref_fasta = repo_dir + '/db/tree_ref.fa'
    args.reference_name = "Wuhan_IPBCAMS-WH-01_2019"
    args.tree = repo_dir + '/db/tree.nwk'
    read1 = os.path.join(args.input_dir, "pipeline", "combined.1.fastq.gz")
    read2 = os.path.join(args.input_dir, "pipeline", "combined.2.fastq.gz")
    subprocess.Popen("minimap2 -ax sr -t 8 %s %s %s | samtools view -b | samtools sort -o %s/var_tree.bam -" %
                     (ref_fasta, read1, read2, working_dir), shell=True).wait()
    subprocess.Popen("samtools mpileup -f %s %s/var_tree.bam > %s" % (ref_fasta, working_dir, args.pileup), shell=True).wait()
profile_dict = get_profile_pileup(args)
variant_dict = get_variants(args)
draw_tree(variant_dict, profile_dict, args.tree, args.reference_name)