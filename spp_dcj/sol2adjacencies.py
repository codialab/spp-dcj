# import from built-in packages
from sys import stdout
import csv

# import from own packages
import spp_dcj.data_utils as du


def cmd_arguments(parser):
    parser.add_argument('sol_file', type=open,
            help='solution file of GUROBI optimizer')
    parser.add_argument('id_to_extremity_map', type=open,
            help='mapping between node IDs and extremities')


def main(args):
    # load data
    # id-to-extremity mapping
    id2ext = dict()
    for line in csv.reader(args.id_to_extremity_map, delimiter = '\t'):
        if line:
            id2ext[line[0]] = tuple(line[1:])

    adjacenciesList, _, weightsDict, _, _, _ = du.parseSOL(args.sol_file, id2ext)

    # write adjacencies
    du.writeAdjacencies(adjacenciesList, weightsDict, stdout)
