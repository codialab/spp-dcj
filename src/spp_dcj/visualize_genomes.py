# import from built-in packages
from sys import stdout

# import from third-party packages
import os#, locale
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')
from matplotlib import pylab as plt
import networkx as nx

# import from own packages
import spp_dcj.data_utils as du

COLOR_MAP = {du.ETYPE_ADJ: 'black', du.ETYPE_ID: 'gray' }
LINESTYLE_MAP = {du.ETYPE_ADJ: 'solid', du.ETYPE_ID: 'dotted' }


def constructGenomeGraph(adjacencies, genes):
    ''' constructs genome graphs from input adjacencies'''

    G = nx.MultiGraph()

    for ext1, ext2 in adjacencies:
        G.add_edge(ext1, ext2, type=du.ETYPE_ADJ)

    for gene in genes:
        # ignore telomeres
        if not gene.startswith('t_'):
            G.add_edge((gene, du.EXTR_HEAD), (gene, du.EXTR_TAIL),
                    type=du.ETYPE_ID)
    return G

def cmd_arguments(parser):
    parser.add_argument('candidateAdjacencies', type=open,
            help='candidate adjacencies of the genomes in the phylogeny')
    parser.add_argument('-i', '--highlight', type=open,
            help='highlight adjacencies in visualization')

def main(args):

    # load data
    candidateAdjacencies = du.parseAdjacencies(args.candidateAdjacencies)

    highlights = None
    if args.highlight:
        highlights = du.parseAdjacencies(args.highlight)['adjacencies']

    # construct genomes and draw figures
    n_species = len(candidateAdjacencies['species'])
    plt.figure(figsize=(10, 5 * (n_species+1)//2))
    for i, (species, adjacencies) in enumerate(candidateAdjacencies[ \
            'adjacencies'].items()):
        # construct genome graph
        genes = candidateAdjacencies['genes'][species]
        G = constructGenomeGraph(adjacencies, genes)
        H = nx.MultiGraph()
        H.add_edges_from(G.edges)

        # plot figure
        plt.subplot((n_species + 1) // 2, 2, i+1)
        layout = nx.spring_layout(G, iterations=200)

        edges = list(G.edges(data = True))

        edge_color, edges = zip(*sorted(zip(list(map(lambda x: \
                COLOR_MAP[x[2]['type']], edges)), edges)))

        if highlights and species in highlights:
            highlight_edges  = list(x for x in highlights[species] if
                    G.has_edge(x[0], x[1]))
            nx.draw_networkx(G, pos=layout, edgelist=highlight_edges, \
                    edge_color = 'red', with_labels=False, node_size=0, \
                    width=4)

        nx.draw_networkx(G, pos=layout, edgelist=edges, edge_color = \
                edge_color, with_labels=False, font_size=1, \
                node_color='orange', node_size=2, width=list(map(lambda x: x \
                == 'black' and 2 or 1, edge_color)))
        plt.title(species)
        ax = plt.gca()
        ax.axis('off')
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

    plt.tight_layout(h_pad=4)
    # plot to output
    with os.fdopen(stdout.fileno(), 'wb', closefd=False) as out:
        plt.savefig(out, format='pdf')
        out.flush()

