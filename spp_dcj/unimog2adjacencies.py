# import from built-in packages
from argparse import ArgumentParser
from sys import stdout

# import from own packages
import spp_dcj.data_utils as du

def cmd_arguments(parser):
    parser.add_argument('unimog_file', type=open,
            help='genomes in UNIMOG format')

def main(args):

    out = stdout
    #
    # load data
    #
    genomes = du.parseUniMoG(args.unimog_file)


    print('\t'.join(('#Species', 'Gene_1', 'Ext_1', 'Species', 'Gene_2',
        'Ext_2', 'Weight')), file = out)
    for genome in genomes:

        gName, _ = genome
        adjs = du.unimog2adjacencies(genome)

        for (g1, ext1), (g2, ext2) in adjs:
            print(f'{gName}\t{g1}\t{ext1}\t{gName}\t{g2}\t{ext2}\t1', file=out)
