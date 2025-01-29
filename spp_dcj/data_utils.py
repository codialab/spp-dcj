# import from built-in packages
from collections import defaultdict, deque
from itertools import product
import logging
import csv
import re
import sys

# import from third-party packages
import networkx as nx

#
# global variables
#

CHR_CIRCULAR    = ')'
CHR_LINEAR      = '|'
EXTR_HEAD       = 'h'
EXTR_TAIL       = 't'
EXTR_CAP        = 'o'

ETYPE_ADJ       = 'adj'
ETYPE_ID        = 'id'
ETYPE_EXTR      = 'extremity'

VTYPE_EXTR      = 'marker_extremity'
VTYPE_CAP       = 'telomere'

ORIENT_NEGATIVE = '-'
ORIENT_POSITIVE = '+'

I_EPSILON=0.05

#separator for ILP variables
SEP='_'

SIGN2EXT_1      = {ORIENT_NEGATIVE:EXTR_TAIL,ORIENT_POSITIVE:EXTR_HEAD}
SIGN2EXT_2      = {ORIENT_NEGATIVE:EXTR_HEAD,ORIENT_POSITIVE:EXTR_TAIL}
EXT2SIGN_1      = {EXTR_TAIL:ORIENT_NEGATIVE,EXTR_HEAD:ORIENT_POSITIVE}
EXT2SIGN_2      = {EXTR_TAIL:ORIENT_POSITIVE,EXTR_HEAD:ORIENT_NEGATIVE}

EXT_COMPLEMENT = {EXTR_HEAD:EXTR_TAIL,EXTR_TAIL:EXTR_HEAD}

TRUE_ADJ_WEIGHT  = 1

DEFAULT_GENE_FAM_SEP = '_'

PAT_ADJ = re.compile(r'^(\w+)@([0-9_]+)$')
PAT_MATCHED_EDGE = re.compile(r'^x(\d+)_(\d+)([^0-9 \t]*) 1\s*$')

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

def complement_id(v):
    (gName, (g, extr)) = v
    return (gName,(g,EXT_COMPLEMENT[extr]))

#
# DATA ACQUISITION, PARSERS
#

def parseTree(data):
    ''' Reads a tree returned as a list of branches [descendant parent] '''
    headerMark = '#'
    delimiter  = '\t'
    res = list()
    for line in csv.reader(data, delimiter = delimiter):
        if line[0][0] != headerMark:
            res.append(line)
    return res


def getTreeDepth(tree):

    # result is a dictionary over all vertices of the tree that maps vertices
    # to their tree depth
    res = dict()

    treeDict = dict(tree)
    revTree = dict()

    root = None
    for child, parent in tree:
        if parent not in revTree:
            revTree[parent] = list()
        revTree[parent].append(child)
        if parent not in treeDict:
            root = parent

    queue = [(0, root)]
    while queue:
        d, v = queue.pop()
        res[v] = d

        if v in revTree:
            for u in revTree[v]:
                queue.append((d+1, u))

    return res

def getFamiliesFromGenes(genesList,speciesList, sep=DEFAULT_GENE_FAM_SEP):
    resFamilies = {}
    for species in speciesList:
        resFamilies[species] = defaultdict(list)
        for gene in genesList[species]:
            family = getFamily(gene, sep=sep)
            resFamilies[species][family].append(gene)
    return(resFamilies)


def addAdjacency(ext1,ext2,weight,resAdjacencies,resWeights,resGenes):
    if ext1>ext2:
        ext1,ext2=ext2,ext1
    gene1 = ext1[0]
    gene2 = ext2[0]
    resAdjacencies.append([ext1,ext2])
    resWeights[ext1,ext2] = weight
    resGenes.append(gene1)
    resGenes.append(gene2)


def parseAdjacencies(data, sep=DEFAULT_GENE_FAM_SEP):
    '''Read a file of adjacencies in the format species\tgene1\text1\tspecies\tgene2\text2\tweight'''
    headerMark = '#'
    delimiter = '\t'
    resAdjacencies = defaultdict(list)
    resGenes       = defaultdict(list)
    resWeights     = defaultdict(dict)
    for line in csv.reader(data, delimiter = delimiter):
        if line[0][0] != headerMark:
            species  = line[0]
            gene1    = line[1]
            ext1     = (gene1,line[2])
            gene2    = line[4]
            ext2     = (gene2,line[5])
            weight   = float(line[6])
            if (gene1, ext1) == (gene2, ext2):
                LOG.warning(f'adjacency between the same extremity is not '
                        f'feasible, skipping {line[1]}:{line[2]}-{line[4]}:{line[5]}')
            else:
                addAdjacency(ext1,ext2,weight,resAdjacencies[species],resWeights[species],resGenes[species])

    speciesList = list(resAdjacencies.keys())
    for species in speciesList:
        resGenes[species]    = list(set(resGenes[species]))
    resFamilies = getFamiliesFromGenes(resGenes, speciesList, sep=sep)

    return {'species':speciesList, 'genes':resGenes,
            'adjacencies':resAdjacencies, 'weights':resWeights,
            'families': resFamilies}


def parseCandidateAdjacencies(data, sep=DEFAULT_GENE_FAM_SEP):
    '''Read candidate adjacencies, returned as a dictionary, indexed by
    species, where for each species we have a list of pairs of gene
    extremities.

    @returns also the species list and the list of genes seen in the
    adjacencies.'''

    headerMark = '#'
    delimiter = ' '
    resAdjacencies = defaultdict(list)
    resGenes       = defaultdict(list)
    resWeights     = defaultdict(dict)
    for line in csv.reader(data, delimiter = delimiter):
        if line[0][0] != headerMark:
            species = line[0]
            m1 = PAT_ADJ.match(line[1])
            m2 = PAT_ADJ.match(line[2])
            if not m1 or not m2:
                raise SyntaxError('unable to parse genes {}/{}'.format(gene1,
                        gene2))
            sp1, gene1   = m1.groups()
            sp2, gene2   = m2.groups()
            if sp1 != sp2 or sp1 != species:
                raise RuntimeError(('adjacencies can only be formed ' + \
                        'between genes of the same species. sp: {} g1: ' + \
                        '{} g2: {}').format(species, line[1], line[2]))
            sign1   = line[3]
            sign2   = line[4]
            weight  = float(line[5])
            ext1 = (gene1,SIGN2EXT_1[sign1])
            ext2 = (gene2,SIGN2EXT_2[sign2])
            addAdjacency(ext1,ext2,weight,resAdjacencies[species],resWeights[species],resGenes[species])
    speciesList = list(resAdjacencies.keys())
    for species in speciesList:
        resGenes[species] = list(set(resGenes[species]))
    resFamilies = getFamiliesFromGenes(resGenes,speciesList, sep=sep)

    return {'species':speciesList, 'genes':resGenes,
            'adjacencies':resAdjacencies, 'weights':resWeights,
            'families': resFamilies}


def parseTrueGeneOrders(data, close_linear=False, sep=DEFAULT_GENE_FAM_SEP):
    '''Read true gene orders, returned as a dictionary, indexed by species,
    where for each species we have a list of pairs of gene extremities.
    If close_linear = True, we assume linear chromosomes and we close them all.
    Assumption: the genes are listed per chromosome and according to their order
    along the chromosome.

    @returns also the species list and the list of genes seen in the
    adjacencies.'''

    headerMark = '#'
    delimiter = '\t'
    resAdjacencies = defaultdict(list)
    resGenes       = defaultdict(list)
    resWeights     = defaultdict(dict)
    prevGene       = ''
    prevSign       = ''
    prevSpecies    = ''
    prevChromosome = ''
    firstGene      = ''
    for line in csv.reader(data, delimiter = delimiter):
        if line[0][0] != headerMark:
            currentSpecies    = line[0]
            currentChromosome = line[1]
            currentGene       = PAT_ADJ.match(line[2]).groups()[1]
            currentSign       = line[3]
            resGenes[currentSpecies].append(currentGene)
            if currentSpecies==prevSpecies and \
              currentChromosome==prevChromosome:
              # we add a new gene to the current chromosome
                ext1 = (prevGene,SIGN2EXT_1[prevSign])
                ext2 = (currentGene,SIGN2EXT_2[currentSign])
                if ext1>ext2:
                    ext1,ext2=ext2,ext1
                resAdjacencies[currentSpecies].append([ext1,ext2])
                resWeights[currentSpecies][ext1,ext2] = TRUE_ADJ_WEIGHT
            else:
                # we start a new chromosome
                if close_linear:
                    if (currentSpecies==prevSpecies and currentChromosome!=prevChromosome) \
                      or \
                      (currentSpecies!=prevSpecies  and prevSpecies!=''):
                      # if close_linear==True and we do not deal with the very first chromosome, we need to close it
                        ext1 = (prevGene,SIGN2EXT_1[prevSign])
                        ext2 = (firstGene,SIGN2EXT_2[firstSign])
                        if ext1>ext2:
                            ext1,ext2=ext2,ext1
                        if [ext1,ext2] not in resAdjacencies[prevSpecies]:
                            resAdjacencies[prevSpecies].append([ext1,ext2])
                            resWeights[prevSpecies][ext1,ext2] = TRUE_ADJ_WEIGHT
                    # We record the current gene as the first gene of the previous chromosome
                    firstGene = currentGene
                    firstSign = currentSign

            prevGene       = currentGene
            prevSign       = currentSign
            prevChromosome = currentChromosome
            prevSpecies    = currentSpecies

    # close last chromosome if anything at all has been read (checked by
    # non-empty prevGene)
    if close_linear and prevGene != '':
        ext1 = (prevGene,SIGN2EXT_1[prevSign])
        ext2 = (firstGene,SIGN2EXT_2[firstSign])
        if ext1>ext2:
            ext1,ext2=ext2,ext1
        if [ext1,ext2] not in resAdjacencies[prevSpecies]:
            resAdjacencies[prevSpecies].append([ext1,ext2])
            resWeights[prevSpecies][ext1,ext2] = TRUE_ADJ_WEIGHT

    speciesList = list(resAdjacencies.keys())
    resFamilies = getFamiliesFromGenes(resGenes,speciesList, sep=sep)

    return {'species':speciesList, 'genes':resGenes,
            'adjacencies':resAdjacencies, 'weights':resWeights,
            'families': resFamilies}


def parseUniMoG(data, genomesOnly=None):
    """Read genome in UniMoG format
    (https://bibiserv.cebitec.uni-bielefeld.de/dcj?id=dcj_manual)"""

    res = list()

    # helper function for parsing each individual gene
    str2gene = lambda x: x.startswith(ORIENT_NEGATIVE) and (ORIENT_NEGATIVE, \
            x[1:]) or (ORIENT_POSITIVE, x.lstrip(ORIENT_POSITIVE))
    # process each line, assuming that the file is well-formatted
    skip = False
    for line in data:
        line = line.strip()
        if line:
            if line.startswith('>'):
                genomeName = line[1:].strip()
                if genomesOnly == None or genomeName in genomesOnly:
                    skip = False
                    res.append((genomeName, list()))
                elif genomesOnly:
                    skip = True
            elif line[-1] not in (CHR_CIRCULAR, CHR_LINEAR):

                raise Exception('Invalid format, expected chromosome to ' + \
                        'end with either \'%s\' or \'%s\'' %(CHR_CIRCULAR, \
                        CHR_LINEAR))
            elif not skip:
                res[-1][1].append((line[-1], list(map(str2gene,
                    line[:-1].split()))))
    return res

def unimog2adjacencies(genome):

    occ = dict()
    res = list()

    # ignore genome name
    _, chromosomes = genome
    tels=0
    for chr_ in chromosomes:
        # set counter for first marker
        if chr_[1][0][1] not in occ:
            occ[chr_[1][0][1]] = 0
        occ[chr_[1][0][1]] += 1

        fst_occ = occ[chr_[1][0][1]]
        for i in range(len(chr_[1])-1):
            (o1, g1), (o2, g2) = chr_[1][i:i+2]
            if g2 not in occ:
                occ[g2] = 0
            res.append(((f'{g1}_{occ[g1]}', SIGN2EXT_1[o1]),
                (f'{g2}_{occ[g2]+1}', SIGN2EXT_2[o2])))
            # increase counter only after-the-fact, in case g1==g2
            occ[g2] += 1

        (o1, g1), (o2, g2) = chr_[1][-1], chr_[1][0]
        if chr_[0] == CHR_CIRCULAR:
            res.append(((f'{g1}_{occ[g1]}', SIGN2EXT_1[o1]),
                (f'{g2}_{fst_occ}', SIGN2EXT_2[o2])))
        elif chr_[0] == CHR_LINEAR:
            tels+=1
            res.append(((f'{g1}_{occ[g1]}', SIGN2EXT_1[o1]), ('t_{}'.format(tels), EXTR_CAP)))
            tels+=1
            res.append((('t_{}'.format(tels), EXTR_CAP), (f'{g2}_{fst_occ}', SIGN2EXT_2[o2])))

    return res

def adjacencies2unimog(adjacenciesList, matchingList):

    genomes = list()
    #
    # assign each family a new unique identifier
    #
    node2fam = lambda x: x[1][:x[1].find('_')]
    famG = nx.Graph(matchingList)
    famC = dict((node2fam(x), 1) for x in famG.nodes() if not
            x[1].startswith('t'))
    for C in nx.connected_components(famG):
        f = node2fam(tuple(C)[0])
        # skip telomeres
        if f not in famC:
            continue
        for v in C:
            famG.nodes[v]['id'] = famC[f]
        famC[f] += 1
    #
    # construct genomes from adjacencies
    for gName, adjs in adjacenciesList.items():
        G = nx.MultiGraph(adjs)
        for adj in adjs:
            for g, ext in adj:
                # iterate through tail extremities, add genes
                if not g.startswith('t') and ext == EXTR_TAIL:
                    G.add_edge((g, EXTR_TAIL), (g, EXTR_HEAD))
        chrs = list()
        for C in nx.connected_components(G):
            degs = set(map(lambda x: x[1], G.degree(C)))
            if degs.issubset((1, 2)):
                isLinear = 1 in degs
                path = None
                if isLinear:
                    v = next((u for u, d in G.degree(C) if d == 1))
                    path = list(nx.traversal.dfs_preorder_nodes(G, v))
                else:
                    v = tuple(C)[0]
                    path = list(nx.traversal.dfs_preorder_nodes(G, v))
                    if len(path) > 2 and path[0][0][0] == path[1][0][0]:
                        path = path[1:] + path[:1]
                chr_ = list()
                for i in range(0, len(path), 2):
                    u = path[i]
                    if u[0].startswith('t'):
                        continue
                    if (gName, u[0]) not in famG:
                        g = f'x_{u[0][:u[0].find("_")]}'
                    else:
                        g = '_'.join((u[0][:u[0].find('_')],
                            str(famG.nodes[(gName, u[0])]['id'])))
                    if u[1] == EXTR_HEAD:
                        chr_.append((ORIENT_POSITIVE, g))
                    elif u[1] == EXTR_TAIL:
                        chr_.append((ORIENT_NEGATIVE, g))
                if chr_:
                    chrs.append((isLinear and CHR_LINEAR or CHR_CIRCULAR, chr_))
                elif not all(map(lambda x: x[0][1:].isdigit(), C)):
                    raise Exception(f'chromosome {C} is empty')
            else:
                raise Exception(f'genome {gName} is not linear/circular')
        genomes.append((gName, chrs))

    return genomes


def parse_edge_id(r):
    return tuple(r.split('_'))

def parseSOL(data, idMap):
    """ SOL file parser """

    obj_value = None
    adjacenciesList = dict()
    matchingList = set()
    matchingDict = dict()
    indelList = dict()
    weightsDict = defaultdict(lambda: defaultdict(float))
    isEmpty = True

    vars_ = dict()

    # objective value is stored in the following comment:
    obj_txt = '# Objective value = '
    for line in data:
        if line.startswith(obj_txt):
            obj_value = float(line[len(obj_txt):])
            continue
        var_, val = line.split()
        vars_[var_] = float(val)
        isEmpty = False
    #collect set adjacencies
    for var_,val in vars_.items():
        if abs(val - 1) > I_EPSILON:
            continue
        if var_.split(SEP)[0]=='a':
            entries = var_.split(SEP)
            #format aSEPedge
            rest = SEP.join(entries[1::])
            id1, id2, _ = parse_edge_id(rest)
            ext1 = idMap[id1]
            ext2 = idMap[id2]
            if ext1[0] not in adjacenciesList:
                adjacenciesList[ext1[0]] = list()
            adj = (ext1[1:], ext2[1:])
            adjacenciesList[ext1[0]].append(adj)
        elif var_.split(SEP)[0]=='x' and not var_.split('_')[-1]=='adj':
            #edge is either indel or match
            entries = var_.split(SEP)
            #format xSEPteSEPedge
            rest = SEP.join(entries[2::])
            id1, id2, etype = parse_edge_id(rest)
            rest = var_[len('a')+len(SEP)::]
            ext1 = idMap[id1]
            ext2 = idMap[id2]
            if etype=='ext':
                e = ext1 < ext2 and (ext1[:2], ext2[:2]) or (ext2[:2],
                            ext1[:2])
                matchingList.add(e)
                matchingDict[ext1] = ext2
                matchingDict[ext1] = ext1
            elif etype=='ind':
                if ext1[0] not in indelList:
                    indelList[ext1[0]] = list()
                indelList[ext1[0]].append((ext1[1:], ext2[1:]))
    return adjacenciesList, indelList, weightsDict, sorted(matchingList), \
            obj_value, vars_

#
# CORE & CONVENIENCE FUNCTIONS
#

def getLeaves(branches):
    '''Creates a boolean dictionary indexed by species where a species has
    value True if it is a leaf'''
    leavesDict = {}
    for [child,parent] in branches:
        leavesDict[parent] = True
        leavesDict[child]  = True
    for [child,parent] in branches:
        leavesDict[parent] = False
    return leavesDict


def getFamily(gene_extr, sep=DEFAULT_GENE_FAM_SEP):
    ''' @returns the family identifier of a gene or gene extremity'''
    assert type(gene_extr) is tuple or type(gene_extr) is str

    # input can either be
    if type(gene_extr) == tuple:
        gene_extr = gene_extr[0]

    return gene_extr[:gene_extr.find(sep)]


def mapFamiliesToGenes(genes, sep=DEFAULT_GENE_FAM_SEP):

    res = dict()
    for gene in genes:
        gid = getFamily(gene, sep=sep)
        if gid not in res:
            res[gid] = list()
        res[gid].append(gene)

    return res


def _constructRDAdjacencyEdges(G, gName, adjacencies, candidateWeights,
        extremityIdManager):
    ''' create adjacencies of the genome named <gName>'''
    for ext1, ext2 in adjacencies:
        id1 = extremityIdManager.getId((gName, ext1))
        id2 = extremityIdManager.getId((gName, ext2))

        # ensure that each edge has a unique identifier
        edge_id = '{}_{}_adj'.format(*sorted((G.nodes[id1]['anc'], G.nodes[id2]['anc'])))
        anc = '{}_{}_adj'.format(*sorted((G.nodes[id1]['anc'],G.nodes[id2]['anc'])))
        weight = candidateWeights.get((ext1, ext2), 0)
        G.add_edge(id1, id2, type=ETYPE_ADJ, id=edge_id, weight=weight,anc=anc)


def _constructNaiveRDCapping(G, gName1, gName2, extremityIdManager):

    caps = dict(((gName1, list()), (gName2, list())))

    for v, vdata in G.nodes(data=True):
        if vdata['type'] == VTYPE_CAP:
            caps[vdata['id'][0]].append(v)

    _addCartesianProductCaps(G, gName1, gName2, caps, extremityIdManager)


def _addCartesianProductCaps(G, gName1, gName2, caps, extremityIdManager):

    new_caps = {gName1: list(), gName2: list()}

    # the number of caps that will be used *at most* is (i) even, and (ii) not
    # more than there are caps in the other genome
    if len(caps[gName1]) < len(caps[gName2]):
        n = (len(caps[gName2])-len(caps[gName1]))//2 * 2
        new_caps[gName1].extend(_fillUpCaps(G, gName1, n, extremityIdManager))
        caps[gName1].extend(new_caps[gName1])
    elif len(caps[gName2]) < len(caps[gName1]):
        n = (len(caps[gName1])-len(caps[gName2]))//2 * 2
        new_caps[gName2].extend(_fillUpCaps(G, gName2, n, extremityIdManager))
        caps[gName2].extend(new_caps[gName2])

    for u, v in product(caps[gName1], caps[gName2]):
        if not G.has_edge(u, v):
            G.add_edge(u, v, type=ETYPE_EXTR, id='{}_{}'.format(*sorted((u, v))))
    return new_caps


def _fillUpCaps(G, gName, ncaps, extremityIdManager):
    new_caps = list()
    for _ in range(ncaps):
        id_ = 't{}'.format(extremityIdManager._IdManager__count)
        v = extremityIdManager.getId((gName, (id_, EXTR_HEAD)))
        new_caps.append(v)
        G.add_node(v, id=(gName, (id_, EXTR_HEAD)), type=VTYPE_CAP)

    for i in range(0, ncaps-1, 2):
        id1 = new_caps[i]
        id2 = new_caps[i+1]
        if not G.has_edge(id1, id2):
            G.add_edge(id1, id2, type=ETYPE_ADJ, id=f'{id1}_{id2}', weight=0.0)
    return new_caps


def _constructRDCapping(G, gName1, gName2, extremityIdManager):

    tel_pairs = _find_end_pairs(G, gName1, gName2)
    A_caps_with_runs, B_caps_with_runs = set(), set()

    # fix paths that are not connected to run-enclosing paths
    norun = {gName1: set(), gName2: set()}
    for u, v, hasArun, hasBrun in tel_pairs:
        if not hasArun and not hasBrun:
            caps = {gName1: list(), gName2: list()}
            caps[G.nodes[u]['id'][0]].append(u)
            caps[G.nodes[v]['id'][0]].append(v)
            _addCartesianProductCaps(G, gName1, gName2, caps,
                    extremityIdManager)
            if G.nodes[u]['id'][0] != G.nodes[v]['id'][0]:
                norun[G.nodes[u]['id'][0]].add(u)
                norun[G.nodes[v]['id'][0]].add(v)
        else:
            for w in (u, v):
                if G.nodes[w]['id'][0] == gName1:
                    A_caps_with_runs.add(w)
                else:
                    B_caps_with_runs.add(w)
    # treat no-run nodes: if their numbers are not matching, add caps to match
    if len(norun[gName1])//2 > len(norun[gName2])//2:
        n = (len(norun[gName1])//2 - len(norun[gName2])//2) * 2
        caps = _fillUpCaps(G, gName2, n, extremityIdManager)
        _addCartesianProductCaps(G, gName1, gName2, \
                {gName1: list(norun[gName1]), gName2: caps}, \
                extremityIdManager)
    elif len(norun[gName2])//2 > len(norun[gName1])//2:
        n = (len(norun[gName2])//2 - len(norun[gName1])//2) * 2
        caps = _fillUpCaps(G, gName1, n, extremityIdManager)
        _addCartesianProductCaps(G, gName1, gName2, \
                {gName1: caps, gName2: list(norun[gName2])}, \
                extremityIdManager)

    _addCartesianProductCaps(G, gName1, gName2, \
            {gName1: list(A_caps_with_runs), gName2: list(B_caps_with_runs)}, \
            extremityIdManager)

#    pos = nx.spring_layout(G)
#    pos = nx.spectral_layout(G)
#
#    genes_edg = list()
#    for gName in {gName1, gName2}:
#        genes = [(u, v) for u, v, data in G.edges(data=True) if data['type'] ==
#                ETYPE_ADJ]
#        genes_edg.extend(((extremityIdManager.getId((gName, (g, EXTR_HEAD))),
#            extremityIdManager.getId((gName, (g, EXTR_TAIL))))for g in genes))
#    Gp = nx.Graph()
#    Gp.add_edges_from(genes_edg)
#    Gp.add_edges_from((u, v) for u, v, data in G.edges(data=True) if
#            data['type'] == ETYPE_ADJ)
#    pos = nx.spring_layout(Gp)
#
#
#    nx.draw_networkx_nodes(G, pos=pos, node_size=2)
#    nx.draw_networkx_labels(G, pos=pos, font_size=8, labels = dict((v,
#        '{0}:{1[0]}{1[1]}'.format(*G.nodes[v]['id'])) for v in G.nodes()))
#
##    nx.draw_networkx_edges(G, pos, [(u, v) for u, v, data in G.edges(data=True)
##        if data['type'] == ETYPE_EXTR], edge_color='green')
#    nx.draw_networkx_edges(G, pos, [(u, v) for u, v, data in G.edges(data=True)
#        if data['type'] == ETYPE_ID], edge_color='gray')
#    nx.draw_networkx_edges(G, pos, [(u, v) for u, v, data in G.edges(data=True)
#        if data['type'] == ETYPE_ADJ], edge_color='red')
#    import matplotlib.pylab as plt
#    plt.savefig('myfig.pdf', format='pdf')

def _find_end_pairs(G, gName1, gName2):
    """ finds all alternating paths between nodes of degree one, which are
    assumed to be caps, i.e., incident to one adjacency edge that connects the
    cap to another node"""

    res = dict()

    # identify caps
    valid_ends = set((v for v, data in G.nodes(data=True) if data['type'] == VTYPE_CAP))

    # checks if edge v-u is ID and if so, sets ID label of corresponding genome
    # to 1, and returns the label vector.
    pID = lambda v, edata, l: G.nodes[v]['id'][0] == gName1 and \
            [l[2] or edata['type'] == ETYPE_ID, l[3]] or \
            [l[2], l[3] or edata['type'] == ETYPE_ID]
    # greater than
    gt = lambda x: x[0] > x[1]

    for end in valid_ends:

        # encoding: state0, state1, has_A_run, has_B_run
        labels = dict(((v, [0, 0, 0, 0]) for v in G.nodes))
        # initialize labeling for root node: caps are connected by
        # adjacency edges (0), so the initial state is 1.
        labels[end][1] = 1
        queue = deque([end])
        while queue:
            v = queue.popleft()
            for u in G.neighbors(v):
                for data in G[u][v].values():
                    # check parity
                    p = data['type'] != ETYPE_ADJ and 1 or 0
                    if labels[v][1-p] > labels[u][p] or (labels[v][1-p] == \
                            labels[u][p] and labels[u][p] and any(map(gt, \
                            zip(pID(v, data, labels[v]), labels[u][2:])))):
                        labels[u][p] = 1
                        labels[u][2] |= labels[v][2]
                        labels[u][3] |= labels[v][3]

                        if G.nodes[u]['id'][0] == gName1:
                            # update A-run flag
                            labels[u][2] |= data['type'] == ETYPE_ID
                        else:
                            # update B-run flag
                            labels[u][3] |= data['type'] == ETYPE_ID

                        if G.nodes[u]['type'] == VTYPE_CAP and u != end:
                            x, y = end < u and (end, u) or (u, end)
                            if (x, y) not in res:
                                res[(x, y)] = [0, 0]
                            res[x, y][0] |= labels[u][2]
                            res[x, y][1] |= labels[u][3]
                        else:
                            queue.append(u)
    return {(u, v, Arun, Brun) for (u, v), (Arun, Brun) in res.items()}


def checkGraph(G,cf=False,checkForAllTels=False):
    for u, v, in G.edges():
        if u == v:
            raise Exception(f'node {v} is connected to itself')

        types = set()
        for data in G[u][v].values():
            if data['type'] not in types:
                types.add(data['type'])
            else:
                raise Exception(f'nodes {u} {G.nodes[u]["id"]}, ' + \
                        f'{v} {G.nodes[v]["id"]} are connected by ' + \
                        f'multiple edges of the type {data["type"]}')

    for v, vdata in G.nodes(data = True):
        hasAdj = False
        hasExtrOrId = False

        if vdata['id'][1][1] not in {EXTR_HEAD, EXTR_TAIL, EXTR_CAP}:
            raise Exception(f'node {v} {G.nodes[v]["id"]} has malformed ' + \
                    'extremity')
        # if vdata['id'][1][1] != EXTR_CAP and checkForAllTels:
        #     has_cap=False
        #     for u in G.neighbors(v):
        #         has_cap=has_cap or G.nodes[u]['type']==EXTR_CAP
        #     if not has_cap:
        #         raise Exception("Looked for caps, but node {} ({}) is lacking one.".format(v,vdata))
        for u in G.neighbors(v):
            for data in G[u][v].values():
                hasAdj |= data['type'] == ETYPE_ADJ
                hasExtrOrId |= data['type'] in {ETYPE_ID, ETYPE_EXTR}
        if not hasAdj:
            raise Exception(f'node {v} {G.nodes[v]["id"]} is not incident ' + \
                    'to an adjacency edge')
        if not hasExtrOrId and (not vdata['type']==VTYPE_CAP or not cf):
            raise Exception(f'node {v} {G.nodes[v]["id"]} is not incident ' + \
                    'to an extremity or indel edge')


class TooManyCSException(Exception):
    pass

def identifyCircularSingletonCandidates(G,max_number=None):
    """ finds all components that can be circular singletons """

    res = dict()
    id_edges = filter(lambda x: x[2]['type'] == ETYPE_ID, G.edges(data=True))
    for e_id in id_edges:
        # orient the traversal: e_id[0] -> e_id[1] -> ...
        # each element of the queue is a tuple of <path, nodeset>
        # - path encoding: <vertex1> <data of edge> <vertex2> ...
        # - node set: set of nodes of the path
        queue = deque((((e_id[0], e_id[2], e_id[1]), set((e_id[0], e_id[1]))), ))
        while queue:
            path, nset = queue.pop()
            v = path[-1]
            # previous edge type
            ptype = path[-2]['type']
            # expected edge type
            etype = ptype == ETYPE_ID and ETYPE_ADJ or ETYPE_ID
            for u in G.neighbors(v):
                for data in G[v][u].values():
                    if data['type'] == etype:
                        if u not in nset:
                            queue.append((path + (data, u), nset.union((u, ))))
                        elif path[0] == u:
                            # no need to check parity, because path is *always*
                            # started with indel edge and no two indel edges
                            # can be adjacent
                            ppath = rotateToMin(path + (data, ))
                            vpath = tuple((ppath[i] for i in range(0,
                                len(ppath), 2)))
                            #TODO: Tell dany I found a bug here!
                            vpath=canonicizePath(vpath)
                            epath = tuple((ppath[i] for i in range(1,
                                len(ppath), 2)))
                            res[vpath] = epath
                            if max_number is not None:
                                if len(res)>max_number:
                                    raise TooManyCSException("Number of CS is over maximum {}.".format(max_number))
    return res


def annotate_v_circ_sing(G):
    #Find all vertices that could be in a circular singleton
    F = G.copy()
    rmedges= [(u,v,k) for u,v,k,etype in G.edges(keys=True,data='type') if etype==ETYPE_EXTR]
    F.remove_edges_from(rmedges)

    for c in nx.connected_components(F):
        S=F.subgraph(c)
        try:
            nx.find_cycle(S)
            for v in c:
                G.nodes[v]['cscandidate']=True
        except nx.NetworkXNoCycle:
            pass


def rotateToMin(path):
    m = min((path[i] for i in range(0, len(path), 2)))
    i = path.index(m)
    return path[i:] + path[:i]

def canonicizePath(path):
    m = min(path)
    i = path.index(m)
    path = path[i:] + path[:i]
    if len(path) > 2:
        if path[1] <= path[-1]:
            return path
        else:
            return (path[0],) + path[-1:0:-1]
    return path




def _constructRDExtremityEdges(G, gName1, gName2, genes, fam2genes1,
        fam2genes2, extremityIdManager):

    genes1 = genes[gName1]
    genes2 = genes[gName2]

    fam2genes = (fam2genes1, fam2genes2)

    fams = set(fam2genes1.keys()).union(fam2genes2.keys())

    #
    # create
    #   - edges between shared extremities of both genomes
    #   - record siblings
    #   - indel edges between families of unequal size
    #
    siblings = list()
    for fam in fams:
        # create extremity edges
        for gene1, gene2 in product(fam2genes1.get(fam, ()), \
                fam2genes2.get(fam, ())):
            id1h = extremityIdManager.getId((gName1, (gene1, EXTR_HEAD)))
            id1t = extremityIdManager.getId((gName1, (gene1, EXTR_TAIL)))
            id2h = extremityIdManager.getId((gName2, (gene2, EXTR_HEAD)))
            id2t = extremityIdManager.getId((gName2, (gene2, EXTR_TAIL)))

            edge_idh = '{}_{}_ext'.format(*sorted((G.nodes[id1h]['anc'], G.nodes[id2h]['anc'])))
            edge_idt = '{}_{}_ext'.format(*sorted((G.nodes[id1t]['anc'], G.nodes[id2t]['anc'])))


            G.add_edge(id1h, id2h, type=ETYPE_EXTR, id=edge_idh)
            G.add_edge(id1t, id2t, type=ETYPE_EXTR, id=edge_idt)

            # ensure sorted order of sibling edges
            siblings.append((edge_idh, edge_idt))

        # create indel edges between genes of smaller family
        for i, gName in enumerate((gName1, gName2)):
            if len(fam2genes[i].get(fam, ())) > len(fam2genes[i-1].get(fam, \
                    ())):
                for gene in fam2genes[i][fam]:
                    idh = extremityIdManager.getId((gName, (gene, EXTR_HEAD)))
                    idt = extremityIdManager.getId((gName, (gene, EXTR_TAIL)))

                    # ensure that each edge has a unique identifier
                    edge_id = '{}_{}_ind'.format(*sorted((G.nodes[idh]['anc'], G.nodes[idt]['anc'])))
                    #if G.has_edge(idh, idt):
                    #    edge_id = '{}_{}'.format(*sorted((idh, idt),
                    #        reverse=True))
                    G.add_edge(idh, idt, type=ETYPE_ID, id=edge_id)

    return siblings

TELOMERE_GENE='t'
def _constructRDNodes(G, gName, genes, globalIdManager,localIdManager):
    ''' create gene extremity nodes for the genome named <gName> '''

    for extr in (EXTR_HEAD, EXTR_TAIL):
        G.add_nodes_from(((localIdManager.getId((gName, (g, extr))),
            dict(id=((gName, (g, extr))), type=VTYPE_EXTR,anc=globalIdManager.getId((gName, (g, extr))))) for g in genes)) #if g!=TELOMERE_GENE))


def _constructRDTelomeres(G, gName, telomeres, extremityIdManager,localIdManager):
    ''' create telomereic extremity nodes for the genome named <gName> '''
    G.add_nodes_from(((localIdManager.getId((gName, (t, 'o'))),
        dict(id=((gName, (t, 'o'))), type=VTYPE_CAP,anc=extremityIdManager.getId((gName, (t, 'o'))))) for t in telomeres))


def hasIncidentAdjacencyEdges(G, v):
    hasAdj = False
    for u in G.neighbors(v):
        for data in G[v][u].values():
            hasAdj = data['type'] == ETYPE_ADJ
            if hasAdj:
                return hasAdj
    return hasAdj


def getIncidentAdjacencyEdges(G, v):
    res = list()
    for u in G.neighbors(v):
        for data in G[v][u].values():
            if data['type'] == ETYPE_ADJ:
                res.append((u, data))
    return res


def constructRelationalDiagrams(tree, candidateAdjacencies, candidateTelomeres,
        candidateWeights, genes, extremityIdManager,
        sep=DEFAULT_GENE_FAM_SEP):
    ''' constructs for each edge of the tree a relational diagram of the
    adjacent genomes'''

    res = dict((('graphs', dict()), ('siblings', dict())))

    for child, parent in tree:
        G = nx.MultiGraph()
        max_tels = len(candidateTelomeres[child]) + len(candidateTelomeres[parent]) +2 #2*sum([len(genes[gnm]) for gnm in [child,parent]])
        localIdManager = IdManager(max_tels)
        for gName in (child, parent):
            _constructRDNodes(G, gName, genes[gName], extremityIdManager,localIdManager)
            _constructRDTelomeres(G, gName, candidateTelomeres[gName],
                                  extremityIdManager,localIdManager)
            _constructRDAdjacencyEdges(G, gName, candidateAdjacencies[gName],
                    candidateWeights[gName], localIdManager)

        fam2genes1 = mapFamiliesToGenes(genes[child], sep=sep)
        fam2genes2 = mapFamiliesToGenes(genes[parent], sep=sep)
        siblings   = _constructRDExtremityEdges(G, child, parent, genes,
                fam2genes1, fam2genes2, localIdManager)

        res['graphs'][(child, parent)] = G
        res['siblings'][(child, parent)] = siblings
        annotate_v_circ_sing(G)


    # create caps at last, assigning them the highest IDs
    #for child, parent in tree:
    #    G = res['graphs'][(child, parent)]
    #    _constructRDCapping(G, child, parent, extremityIdManager)
#        _constructNaiveRDCapping(G, child, parent, extremityIdManager)
        # remove caps from the graph that are not saturated by two edges
#        for v, d in tuple(G.degree()):
#            if d == 1:
#                if G.nodes[v]['type'] != VTYPE_CAP:
#                    raise Exception('assumption that all vertices with ' + \
#                            'degree one are caps failed')
#                G.remove_node(v)

    return res


def writeAdjacencies(adjacenciesList, weightsDict, out):
    ''' Write an adjacency file '''
    out.write('#Species Gene_1 Ext_1 Species Gene_2 Ext_2 Weight\n')
    speciesList = adjacenciesList.keys()
    for species in speciesList:
        for [(gene1,ext1),(gene2,ext2)] in adjacenciesList[species]:
            out.write('\t'.join([
                species,
                gene1,
                ext1,
                species,
                gene2,
                ext2,
                str(weightsDict[species][(gene1,ext1), (gene2,ext2)])])+'\n')


#
# DATA CLASSES
#

class IdManager(object):

    def __init__(self,cap_limit,is_cap = lambda x : x[1][1]=='o'):
        self.__count_caps = 1
        self.__count_ext = cap_limit+1
        self.__cap_limit = cap_limit
        self.__table = dict()
        self.__reverse = dict()
        self.__is_cap = is_cap

    def getId(self, obj):
        if obj not in self.__table:
            if self.__is_cap(obj):
                if self.__count_caps + 1 > self.__cap_limit:
                    raise AssertionError("Caps have overflowed.")
                self.__table[obj] = self.__count_caps
                self.__reverse[self.__count_caps] = obj
                self.__count_caps += 1
            else:
                self.__table[obj] = self.__count_ext
                self.__reverse[self.__count_ext] = obj
                self.__count_ext += 1
        return self.__table[obj]

    def getObj(self, id_):
        return self.__reverse[id_-1]

    def getMap(self):
        return dict(self.__table.items())
