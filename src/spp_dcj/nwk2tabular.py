# import from built-in packages
from sys import stdout

# import from own packages
from spp_dcj.newick_parser import parse_tree


def cmd_arguments(parser):
    parser.add_argument('tree', type=open, help='tree in newick format')


def main(args):

    # load data
    tree = parse_tree(args.tree)

    out = stdout
    print('#Node\tParent', file = out)
    for v in tree.getNodes():
        u = v.getAncestor()
        if u != None:
            print('\t'.join((v.label, u.label)), file = out)


