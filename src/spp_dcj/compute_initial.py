# import from built-in packages
from sys import stdout

# import from own packages
import spp_dcj.data_utils as du

def cmd_arguments(parser):
    parser.add_argument('tree', type=open,
            help='phylogenetic tree as child->parent relation table (Important: children must be in first column.)')
    parser.add_argument('candidateAdjacencies', type=open,
            help='candidate adjacencies of the genomes in the phylogeny')
    parser.add_argument('idmap',help="Id map file (output of the main script)",type=open)
    parser.add_argument('-t', '--no_telomeres', action='store_true',
            help='don\'t add any additional telomeres. Overrides -at.')
    parser.add_argument('-l', '--all_telomeres', action='store_true',
            help='Add additional telomeres to all nodes.')
    parser.add_argument('--family-bounds','-fmb',type=open)
    parser.add_argument('--set-circ-sing-handling',choices=['adaptive','enumerate','barbershop'],default='adaptive')
    parser.add_argument('-s', '--separator', default = du.DEFAULT_GENE_FAM_SEP, \
            help='Separator of in gene names to split <family ID> and ' +
                    '<uniquifying identifier> in adjacencies file')
    parser.add_argument('-a', '--alpha', default = 0.5, type=float,
            help='linear weighting factor for adjacency weights vs DCJ ' + \
                    'indel distance (alpha = 1 => maximize only DCJ indel ' + \
                    'distance)')
    parser.add_argument('--def-telomere-weight', '-w', default = 0, type=float,
            help='Default weight for added telomeres. Has no effect if -t is used. For most use cases this should be <= 0.')
    parser.add_argument('-b', '--beta', type=float,
            help='Backwards compatible beta parameter from v1. Telomeres will be re-weighted and the ILP scaled to be equivalent to v1.')


def main(args):
    out = stdout

    # load data
    glob_table, loc_tables = du.parse_id_file(args.idmap)
    global_ext2id = du.IdManager(0,is_cap=lambda x: False,initial_table=glob_table)
    speciesTree = du.parseTree(args.tree)

    candidateAdjacencies = du.parseAdjacencies(args.candidateAdjacencies,
                                               sep=args.separator)
    leaves = set([x for x, v in du.getLeaves(speciesTree).items() if v])
    # add telomeres
    telomeres = du.identifyCandidateTelomeres(candidateAdjacencies,args.def_telomere_weight,leaves, dont_add=args.no_telomeres,addToAll=args.all_telomeres)

    # construct adjacency graphs
    genes = candidateAdjacencies['genes']
    adjacencies = candidateAdjacencies['adjacencies']
    weights = candidateAdjacencies['weights']
    fam_bounds = dict()
    if args.family_bounds:
        fam_bounds=du.parseFamilyBounds(args.family_bounds)
    families = candidateAdjacencies['families']
    du.fillFamilyBounds(families,fam_bounds)
    relationalDiagrams = du.constructRelationalDiagrams(speciesTree,
            adjacencies, telomeres, weights, genes, global_ext2id,fam_bounds=fam_bounds,
            sep=args.separator,loc_manager_tables=loc_tables)
    graphs = relationalDiagrams['graphs']
    siblings = relationalDiagrams['siblings']
    for G in graphs.values():
        du.checkGraph(G,cf=True,checkForAllTels=args.all_telomeres and not args.no_telomeres)
    circ_singletons = dict()
    for ident, G in graphs.items():
        if args.set_circ_sing_handling=='adaptive':
            max_tolerable = len([u for u in G.nodes() if G.nodes[u].get('cscandidate',False)])
        elif args.set_circ_sing_handling=='enumerate':
            max_tolerable = None
        else:
            max_tolerable=0
        if max_tolerable==None:
            circ_singletons[ident] = du.identifyCircularSingletonCandidates(G)
        elif max_tolerable > 0:
            try:
                circ_singletons[ident] = du.identifyCircularSingletonCandidates(G,max_number=max_tolerable)
            except du.TooManyCSException:
                pass
    if args.beta:
        scale_factor = 1-args.beta
        our_alpha = args.alpha/scale_factor
        our_beta_add_telweight=args.beta/(scale_factor*(our_alpha-1))
    else:
        scale_factor=1
        our_alpha=args.alpha
        our_beta_add_telweight=0

    if args.beta:
        #rescale the edge weights for the telomeres
        for ident,G in graphs.items():
            du.add_weight_tels(G,our_beta_add_telweight)
    du.warm_start_decomposition(graphs,fam_bounds,families,siblings)
    du.sol_from_decomposition(graphs,circ_singletons,our_alpha, out)
