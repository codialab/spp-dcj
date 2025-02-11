# builtin imports
from argparse import ArgumentParser
import logging
import sys

# own module imnports
from spp_dcj import \
        compute_initial as heuristic, \
        construct_ilp as ilp, \
        sol2adjacencies as sol2adj, \
        nwk2tabular as newick2table, \
        visualize_genomes as draw, \
        __version__ as VERSION

TOOLNAME = 'spp-dcj'

st = logging.StreamHandler(sys.stderr)
LOG = logging.getLogger(TOOLNAME)
LOG.setLevel(logging.DEBUG)

st.setLevel(logging.DEBUG)
st.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
LOG.addHandler(st)

def main():
    parser = ArgumentParser(prog=TOOLNAME,
                            description='A tool to compute the a solution to the small parsimon problem of natural genomes')
    parser.add_argument('-v', '--version', help='show version and exit', action='version', version=f'{TOOLNAME} {VERSION}')

    subparsers = parser.add_subparsers(dest='subparser')
    heuristicparser = subparsers.add_parser('heuristic',
                                      help='computes heuristic solution for SPP-DCJ in SOL format that can be used as warm-start for the ILP')
    heuristic.cmd_arguments(heuristicparser)
    ilpparser = subparsers.add_parser('ilp',
                                      help='construct an ILP for computing a solution to the small parsimony problem of natural genomes')
    ilp.cmd_arguments(ilpparser)
    sol2adjparser = subparsers.add_parser('sol2adj',
                                      help='converts a GUROBI solution file to a tab-separated table of adjacencies')
    sol2adj.cmd_arguments(sol2adjparser)
#    sol2uniparser = subparsers.add_parser('sol2unimog',
#                                      help='converts a GUROBI solution file to a genome file in UniMoG format')
#    sol2unimog.cmd_arguments(sol2uniparser)

    nwk2tabparser = subparsers.add_parser('newick2table',
                                      help='translate tree in Newick format to tabular format that is required by SPP-DCJ ILP')
    newick2table.cmd_arguments(nwk2tabparser)
#    unimog2adjparser = subparsers.add_parser('unimog2adj',
#                                      help='translate genomes in UniMoG format to adjacencies list requires as input to SPP-DCJ ILP')
#    unimog2adj.cmd_arguments(unimog2adjparser)
    visualparser = subparsers.add_parser('draw',
                                      help='visualize one or more degenerate genomes as graph')
    draw.cmd_arguments(visualparser)
    if len(sys.argv) < 2:
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    globals()[args.subparser].main(args)

