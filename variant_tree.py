import argparse, sys
from ete3 import Tree, RectFace, AttrFace, TextFace
from ete3 import NodeStyle
from ete3 import TreeStyle

def get_profile_pileup(args):
    outlist = []
    with open(args.pileup) as f:
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
            if depth >= args.min_depth:
                variants = []
                for i in ['a', 't', 'c', 'g']:
                    if counts[i] /depth > args.min_fraction:
                        variants.append((depth, i))
                if len(variants) >= 2:
                    variants.sort(reverse=True)
                    outlist.append([pos, []])
                    for i in variants:
                        outlist[1].append(i[1])
                elif len(variants) == 1 and variants[0][1] != refbase and args.include_variants:
                    outlist.append([pos, [variants[0][1]]])
    return outlist





def get_array(args, variant_sites):
    outdict = {}
    seq_dict = {}
    ref_name = args.reference_name
    with open(args.multiple_alignment) as f:
        for line in f:
            if line.startswith(">"):
                name = line.split()[0][1:]
                seq_dict[name] = ''
            else:
                seq_dict[name] += line.rstrip().lower()
    pos_translate = {}
    count = 1
    for i in range(0, len(seq_dict[ref_name])):
        if seq_dict[ref_name] == '-':
            pass
        else:
            pos_translate[i] = count
            count += 1

    for i in variant_sites:
        pos = pos_translate[i[0]]
        for j in seq_dict:
            if not j in outdict:
                outdict[j] = {}
            if seq_dict[j][pos] == variant_sites[1][0]:
                outdict[j][pos] = 1
            elif seq_dict[j][pos] in variant_sites[1]:
                outdict[j][pos] = 2
            else:
                outdict[j][pos] = 0
    return outdict



    variant_dict = {}
    return variant_dict

def draw_tree(variant_dict, the_tree):
    t = Tree(the_tree, quoted_node_names=True)
    font_size = 8
    font_type = 'Heveltica'
    font_gap = 3
    font_buffer = 10
    t.set_outgroup()
    ts = TreeStyle()
    ts.legend.add_face(
        TextFace(font_gap * ' ' + 'Gene present/absent' + ' ' * font_buffer, fsize=font_size, ftype=font_type,
                 tight_text=True), column=starting_col + 1)
    ts.legend.add_face(RectFace(20, 20, '#FFFFFF', '#FFFFFF'), column=0)
    ts.legend.add_face(
        TextFace(font_gap * ' ' + 'Gene present/absent' + ' ' * font_buffer, fsize=font_size, ftype=font_type,
                 tight_text=True), column=1)
    ts.legend.add_face(RectFace(20, 20, "#5ba965", "#5ba965"), column=0)
    ts.legend.add_face(
        TextFace(font_gap * ' ' + 'Gene present/absent' + ' ' * font_buffer, fsize=font_size, ftype=font_type,
                 tight_text=True), column=1)
    ts.legend.add_face(RectFace(20, 20, "#cb5b4c", "#cb5b4c"), column=0)

    nameF = TextFace(font_gap * ' ' + name + ' ' * font_buffer, fsize=font_size, ftype=font_type, tight_text=True)
    nameF.rotation = -90
    ts.aligned_header.add_face(nameF, column=starting_col + len(gene_list) - 1)
    for n in t.iter_leaves():
        the_name = n.name
        if the_name[0] == '"' and the_name[-1] == '"':
            the_name = the_name[1:-1]
        if the_name.endswith('.ref'):
            the_name = the_name[:-4]
        for num, i in enumerate(gene_list):
            if i in gotit:
                colour = "#5ba965"
            else:
                colour = "#cb5b4c"
            n.add_face(RectFace(20, 20, colour, colour), column=num + starting_col, position="aligned")


    t.show(tree_style=ts)



__version__ = "0.0.1"





parser = argparse.ArgumentParser(prog='Variant_tree.py', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='Show where hetrozygous bases show up on a tree.\n' \
                                            'Version: %s\n' 
                                            'License: GPLv3\n'
                                            'USAGE: python -p output.pileup -m multiple_alignment -r reference_name -t tree.nwk' % __version__)


parser.add_argument('-m', '--multiple_alignment', action='store', default="12", help='number of threads to use')
parser.add_argument('-r', '--reference_name', action='store', help='number of threads to use')
parser.add_argument('-t', '--tree', action='store', help='number of threads to use')
parser.add_argument('-p', '--pileup', action='store', help='pileup of reads')
parser.add_argument('-f', '--min_fraction', action='store', default=0.05, help='minimum fraction of alternate to flag')
parser.add_argument('-d', '--min_depth', action='store', default=40, help='minimum depth to consider site')
parser.add_argument('-v', '--include_variants', action='store_true')



args = parser.parse_args()
variant_sites = get_profile_pileup(args)
variant_dict = get_array(args, variant_sites)
print(variant_dict)
draw_tree(variant_dict, args.tree)