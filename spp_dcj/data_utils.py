# import from built-in packages
from collections import defaultdict, deque
from itertools import product
from functools import reduce
from math import ceil
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

def get_genome(G, v):
    return G.nodes[v]['id'][0]

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
    '''Read genome in UniMoG format
    (https://bibiserv.cebitec.uni-bielefeld.de/dcj?id=dcj_manual)'''

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
                        g = f'x_{u[0][:u[0].find('_')]}'
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
    ''' SOL file parser '''
    raise NotImplementedError('Not implemented for this version of SPP-DCJ')
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

def parseSOLAdj(data, idMap):
    ''' SOL file parser '''
    obj_value = None
    adjacenciesList = dict()


    vars_ = dict()

    # objective value is stored in the following comment:
    obj_txt = '# Objective value = '
    for line in data:
        if line.startswith(obj_txt):
            obj_value = float(line[len(obj_txt):])
            continue
        #print(line,file=sys.stderr)
        var_, val = line.split()
        vars_[var_] = float(val)
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
            adj = (ext1[1], ext2[1])
            adjacenciesList[ext1[0]].append(adj)
    return adjacenciesList, obj_value


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


def nodeGetFamily(G,v):
    return getFamily(G.nodes[v]['id'][1][0])

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
    ''' finds all alternating paths between nodes of degree one, which are
    assumed to be caps, i.e., incident to one adjacency edge that connects the
    cap to another node'''

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
    #for v,data in G.nodes(data=True):
    #    print(v,','.join(['{}={}'.format(k,x) for k,x in data.items()]),file=sys.stderr)
    gnm_min = dict()
    gnm_max = dict()
    for v in G.nodes():
        gnm = get_genome(G,v)
        if not gnm in gnm_min:
            gnm_min[gnm]=v
            gnm_max[gnm]=v
        gnm_min[gnm]= min(gnm_min[gnm],v)
        gnm_max[gnm]= max(gnm_min[gnm],v)
    gnms = list(gnm_min.keys())
    assert(len(gnms)==2)
    assert(gnm_min[gnms[0]]>gnm_max[gnms[1]] or gnm_min[gnms[1]]>gnm_max[gnms[0]])
    for u, v, in G.edges():
        if u == v:
            raise Exception(f'node {v} is connected to itself')

        types = set()
        for data in G[u][v].values():
            if data['type'] not in types:
                types.add(data['type'])
            else:
                raise Exception(f'nodes {u} {G.nodes[u]['id']}, ' + \
                        f'{v} {G.nodes[v]['id']} are connected by ' + \
                        f'multiple edges of the type {data['type']}')

    for v, vdata in G.nodes(data = True):
        hasAdj = False
        hasExtrOrId = False

        if vdata['id'][1][1] not in {EXTR_HEAD, EXTR_TAIL, EXTR_CAP}:
            raise Exception(f'node {v} {G.nodes[v]['id']} has malformed ' + \
                    'extremity')
        # if vdata['id'][1][1] != EXTR_CAP and checkForAllTels:
        #     has_cap=False
        #     for u in G.neighbors(v):
        #         has_cap=has_cap or G.nodes[u]['type']==EXTR_CAP
        #     if not has_cap:
        #         raise Exception('Looked for caps, but node {} ({}) is lacking one.'.format(v,vdata))
        for u in G.neighbors(v):
            for data in G[u][v].values():
                hasAdj |= data['type'] == ETYPE_ADJ
                hasExtrOrId |= data['type'] in {ETYPE_ID, ETYPE_EXTR}
        if not hasAdj:
            raise Exception(f'node {v} {G.nodes[v]['id']} is not incident ' + \
                    'to an adjacency edge')
        if not hasExtrOrId and (not vdata['type']==VTYPE_CAP or not cf):
            raise Exception(f'node {v} {G.nodes[v]['id']} is not incident ' + \
                    'to an extremity or indel edge')


class TooManyCSException(Exception):
    pass

def identifyCircularSingletonCandidates(G,max_number=None):
    ''' finds all components that can be circular singletons '''

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
                                    raise TooManyCSException('Number of CS is over maximum {}.'.format(max_number))
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
        fam2genes2, extremityIdManager,fam_bounds):

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

        if not fam in fam_bounds[gName1] and not fam in fam_bounds[gName2]:
            raise AssertionError('WTF')
        # create indel edges between genes of smaller? (leobkmer says: do you mean larger?) family
        for i, gName in enumerate((gName1, gName2)):
            oName = [oName for oName in (gName1,gName2) if oName!=gName].pop()
            if fam_bounds[gName].get(fam,(0,0))[1] > fam_bounds[oName].get(fam,(0,0))[0]:
                #print('Fam {} in genome {} overrepresented. Adding indel edges...'.format(fam,gName),file=sys.stderr)
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
        candidateWeights, genes, extremityIdManager,fam_bounds=dict(),
        sep=DEFAULT_GENE_FAM_SEP,loc_manager_tables=dict()):
    ''' constructs for each edge of the tree a relational diagram of the
    adjacent genomes'''

    res = dict((('graphs', dict()), ('siblings', dict()),('idmanagers',dict())))
    for child, parent in tree:
        G = nx.MultiGraph()
        max_tels = len(candidateTelomeres[child]) + len(candidateTelomeres[parent]) +2 #2*sum([len(genes[gnm]) for gnm in [child,parent]])
        if (child,parent) in loc_manager_tables:
            localIdManager=IdManager(max_tels,initial_table=loc_manager_tables[(child,parent)])
        else:
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
                fam2genes1, fam2genes2, localIdManager,fam_bounds)
        
        res['graphs'][(child, parent)] = G
        res['siblings'][(child, parent)] = siblings
        res['idmanagers'][(child,parent)]=localIdManager
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


def getNotFullFamilies(fam_bounds,families):
    not_full_fams = dict()
    for g in fam_bounds:
        if not g in not_full_fams:
            not_full_fams[g] = dict()
        for f in fam_bounds[g]:
            if len(families[g][f]) > fam_bounds[g][f][1]:
                not_full_fams[g][f]=len(families[g][f]) - fam_bounds[g][f][1]
    return not_full_fams

    

#
# DATA CLASSES
#

class IdManager(object):

    def __init__(self,cap_limit,is_cap = lambda x : x[1][1]=='o',initial_table=None):
        self.__count_caps = 1
        self.__count_ext = cap_limit+1
        self.__cap_limit = cap_limit
        self.__table = dict()
        self.__reverse = dict()
        self.__is_cap = is_cap
        if initial_table is not None:
            for k,v in initial_table.items():
                self.__table[k]=v
                self.__reverse[v]=k
    def getId(self, obj):
        if obj not in self.__table:
            if self.__is_cap(obj):
                if self.__count_caps + 1 > self.__cap_limit:
                    raise AssertionError('Caps have overflowed.')
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


def parse_lower_bound_file(lbf):
    lb_map = {}
    with open(lbf) as f:
        for line in f:
            entries = line.strip().split()
            if len(entries)==0:
                continue
            assert(len(entries)==3)
            a,b,lb=entries
            if not a in lb_map:
                lb_map[a]=dict()
            try:
                lb_map[a][b]=float(lb)
            except ValueError:
                pass
    return lb_map


def create_adjacency_map(adjacency_list):
    a_map=dict()
    for x,y in adjacency_list:
        assert(x not in a_map)
        assert(y not in a_map)
        a_map[x]=y
        a_map[y]=x
    return a_map


def invert_extremity(yx):
    y,x = yx
    x_ = EXTR_HEAD if x==EXTR_TAIL else EXTR_TAIL
    return y,x_

def get_family_map(families):
    f_map = dict()
    for f,gs in families.items():
        for g in gs:
            assert(g not in f_map)
        f_map[g]=f
    return f_map

def get_unimog(genome,genes,adjacencies,sep):
    a_map=create_adjacency_map(adjacencies[genome])
    telomeres = []
    visited = dict()
    for gene in genes[genome]:
        for xt in [EXTR_HEAD,EXTR_TAIL]:
            extremity = (gene,xt)
            visited[extremity]=False
            if extremity not in a_map:
                telomeres.append(extremity)
    assert(len(telomeres)%2==0)
    lin_chroms = []
    for t in telomeres:
        curr_chrom=[]
        if visited[t]:
            continue
        curr = t
        
        while True:
            assert(not visited[curr])
            visited[curr]=True
            gene,xt = invert_extremity(curr)
            assert(not visited[(gene,xt)])
            visited[(gene,xt)]=True
            if xt==EXTR_HEAD:
                curr_chrom.append(getFamily(gene,sep))
            else:
                curr_chrom.append('-'+getFamily(gene,sep))
            if (gene,xt) in a_map:
                curr = a_map[(gene,xt)]
            else:
                break
        lin_chroms.append(curr_chrom)
    circ_chroms = []
    for xtr in visited:
        if visited[xtr]:
            continue
        curr = xtr
        curr_chrom=[]
        while True:
            assert(not visited[curr])
            visited[curr]=True
            gene,xt = invert_extremity(curr)
            assert(not visited[(gene,xt)])
            visited[(gene,xt)]=True
            if xt==EXTR_HEAD:
                curr_chrom.append(getFamily(gene,sep))
            else:
                curr_chrom.append('-'+getFamily(gene,sep))
            assert((gene,xt) in a_map)
            curr = a_map[(gene,xt)]
            if visited[curr]:
                break
        circ_chroms.append(curr_chrom)
    unimog_str = '>{}'.format(genome)
    for l in lin_chroms:
        unimog_str+='\n'
        unimog_str+=(' '.join(l))
        unimog_str+=' |'
    for c in circ_chroms:
        unimog_str+='\n'
        unimog_str+=(' '.join(c))
        unimog_str+=' )'
    return unimog_str

def parseFamilyBounds(data):
    bounds = dict()
    delimiter = '\t'
    for line in csv.reader(data, delimiter = delimiter):
        genome,fam,low,high = line
        if genome not in bounds:
            bounds[genome]=dict()
        if fam in bounds[genome]:
            LOG.warning('bound for family {} set twice for genome {}, will be overwritten.'.format(fam,genome))
        bounds[genome][fam]=(int(low),int(high))
    return bounds

def fillFamilyBounds(families,bounds):
    for genome in families:
        if not genome in bounds:
            bounds[genome]=dict()
        for f in families[genome]:
            if not f in bounds[genome]:
                #set both to maximum
                fsize = len(families[genome][f])
                bounds[genome][f]=(fsize,fsize)


def parse_id_file(id_file):
    local_tables = dict()
    global_table = dict()
    for line in csv.reader(id_file, delimiter = '\t'):
        assert(len(line) in [4,6])
        if len(line)==4:
            #global id
            #(str(v), k[0], k[1][0], k[1][1])
            obj = (line[1],(line[2],line[3]))
            oid = int(line[0])
            global_table[obj]=oid
        else:
            pair = (line[0],line[1])
            if not pair in local_tables:
                local_tables[pair]=dict()
            oid = int(line[2])
            obj = (line[3],(line[4],line[5]))
            local_tables[pair][obj]=oid
    return (global_table,local_tables)


def identifyCandidateTelomeres(candidateAdjacencies, telomere_default_weight,leaves, dont_add=False,addToAll=False):

    res = dict()
    weights = candidateAdjacencies['weights']
    for species, adjs in candidateAdjacencies['adjacencies'].items():
        genes = candidateAdjacencies['genes'][species]
        # add gene extremities incident to telomeric adjacencies to telomere
        # set
        telomeres = set(x[0][0] for x in adjs if x[0][1] == 'o').union(
                (x[1][0] for x in adjs if x[1][1] == 'o'))

        # remove telomeric extremities from gene set
        for t in telomeres:
            genes.remove(t)

        if not dont_add:
            G = nx.Graph()
            G.add_nodes_from(reduce(lambda x, y: x + y, (((g, EXTR_HEAD), (g,
                EXTR_TAIL)) for g in genes)))
            G.add_edges_from(adjs)

            for C in tuple(nx.connected_components(G)):
                C = set(C)

                # check if component is linear / circular / fully connected /
                # even odd
                # - if it is (linear/circular or fully connected) and even, no
                # telomere needs to be added
                degs = set(map(lambda x: x[1], G.degree(C)))

                # structures that support a perfect matching, e.g.
                # - simple paths/simple cycles of even size
                # - fully connected components of even size
                # do not need telomeres and are omitted by the condition below.
                # Conversely, nodes of components of odd size (including components
                # of size 1) will always be considered as candidate telomeres

                # will evaluate to true if component is NOT
                # - linear/circular or
                # - fully connected
                # - even (in terms of #nodes)
                if degs.difference((1, 2)) and degs.difference((len(C)-1,)) or len(C) % 2 or (addToAll and species not in leaves):
                    for g, extr in C:
                        #skip already existing telomeres
                        if extr=='o':
                            continue
                        t = f't_{g}_{extr}'
                        telomeres.add(t)
                        adjs.append(((g, extr), (t, 'o')))
                        candidateAdjacencies['weights'][((g, extr), (t, 'o'))]=telomere_default_weight

#        genes_edg = [((g, EXTR_HEAD), (g, EXTR_TAIL)) for g in genes]
#        if species == 'n3':
#            C = nx.connected.node_connected_component(G, ('69_6', 'h'))
#            G = G.subgraph(C).copy()
##            G.add_edges_from(genes_edg)
#            pos = nx.spring_layout(G)
#            #nx.draw_networkx_nodes(G, pos=pos, node_size=8)
#            nx.draw(G, pos=pos, node_size=10)
#            #nx.draw_networkx_labels(G, pos=pos, font_size=12)
#            nx.draw_networkx_edges(G, pos, set(map(tuple,
#                adjs)).intersection(G.edges()), edge_color='black')
##            nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=dict(((x[0], \
##                    x[1]), G[x[0]][x[1]][0]['type']) for x in G.edges(data = \
##                    True)))
##            nx.draw_networkx_edges(G, pos, set(genes_edg).intersection(G.edges()),
##                edge_color='red')
#            import matplotlib.pylab as plt
#            import pdb; pdb.set_trace()
        res[species] = telomeres
    return res


def greedy_matching_completion(G,sibmap):
    for v in G.nodes():
        extren = [(v,u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==ETYPE_EXTR]
        iden = [(v,u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==ETYPE_ID]
        set_extren = [(x,y,z) for x,y,z in extren if G[x][y][z].get('is_set',False)]
        set_iden = [(x,y,z) for x,y,z in iden if G[x][y][z].get('is_set',False)]
        assert(len(set_extren)<=1)
        assert(len(set_iden)<=1)
        assert(len(set_iden)+len(set_extren)<=1)
        if G.nodes[v]['type']==VTYPE_CAP:
            #skip caps, they dont need matching
            continue
        if len(set_extren)==1 or len(set_iden)==1:
            continue
        possibilities = [(x,y,z) for x,y,z in extren if G[x][y][z].get('is_set',True)]
        if len(possibilities)==0:
            #needs deletion
            assert(len(iden)==1)
            x,y,z = iden.pop()
            G[x][y][z]['is_set']=True
        else:
            v,u,k = possibilities.pop()
            #get the sibling
            v_,u_ = sibmap[(v,u)]
            k__candidates = [k_ for k_ in G[v_][u_] if G[u_][v_][k_]['type']==ETYPE_EXTR]
            assert(len(k__candidates)==1)
            k_ = k__candidates.pop()
            corr_at = {u:v,v:u,v_:u_,u_:v_}
            wrongs = set()
            for x,y in corr_at.items():
                wrong = set([(min(w,x),max(w,x),k) for w in G[x] for k in G[x][w] if G[x][w][k]['type']==ETYPE_EXTR and w!=y])
                wrongs.update(wrong)
            for x,y,z in wrongs:
                assert(not G[x][y][z].get('is_set',False))
                G[x][y][z]['is_set']=False
            G[u][v][k]['is_set']=True
            G[u_][v_][k_]['is_set']=True


def check_connectedness(G,v_):
    #check the family of a node and determine if it obeys MMM and is in the smaller family
    gnma = {v_}
    gnmb = set()
    changed = True
    while changed:
        changed = False
        for v in gnma:
            for u in G[v]:
                for k in G[v][u]:
                    if G[v][u][k]['type']==ETYPE_EXTR and u not in gnmb:
                        gnmb.add(u)
                        changed=True
        for v in gnmb:
            for u in G[v]:
                for k in G[v][u]:
                    if G[v][u][k]['type']==ETYPE_EXTR and u not in gnma:
                        gnma.add(u)
                        changed=True
    assert(len(gnma)<=len(gnmb))
    for x in gnma:
        for y in gnmb:
            xte_candidates = [k for k in G[x][y] if G[x][y][k]['type']==ETYPE_EXTR]
            assert(len(xte_candidates)==1)


def fill_matching_greedy(G,sibmap):
    for v in G.nodes():
        exten = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==ETYPE_EXTR and G[v][u][k].get('is_set',False)]
        assert(len(exten)<=1)
        adjen = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==ETYPE_ADJ and G[v][u][k].get('is_set',False)]

        iden = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==ETYPE_ID and G[v][u][k].get('is_set',False)]

        assert(len(adjen)==1)
        assert(len(iden)<=1)
        if G.nodes[v]['type']==VTYPE_CAP and len(adjen) == 0:
            #non-active telomere, skip
            continue
        if len(exten) == 0 and len(iden)==0:
            #not matched yet
            if G.nodes[v]['type']==VTYPE_CAP:
                continue
            success=False
            for u in G[v]:
                for k in G[v][u]:
                    xt_cands_u = [(x,k) for x in G[u] for k in G[u][x] if G[u][x][k]['type']==ETYPE_EXTR and G[u][x][k].get('is_set',False)]
                    if G[v][u][k]['type']==ETYPE_EXTR and G[v][u][k].get('is_set',True) and len(xt_cands_u)==0:
                        G[v][u][k]['is_set']=True
                        v_,u_ = sibmap[(v,u)]
                        candidates = [k for k in G[v_][u_]]
                        #extremity edges are always unique
                        assert(len(candidates)==1)
                        G[v_][u_][candidates[0]]['is_set']=True
                        success=True
                        break
            if not success:
                iden = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==ETYPE_ID]
                try:
                    assert(len(iden)==1)
                except AssertionError:
                    check_connectedness(G,v)
                    LOG.error('connectedness check failed after greedy matching')
                    raise AssertionError('connectedness check failed after greedy matching')
                #set the indel edge
                u,k=iden[0]
                G[v][u][k]['is_set']=True


def set_cycles(G,sibmap,cycles):
    for cyc in cycles:
        to_set = set()
        to_remove = set()
        for u,v,k in cyc:
            if v < u:
                u,v=v,u
            u_,v_ = sibmap[(u,v)]
            if v_ < u_:
                u_,v_ = v_,u_
            cand = [k_ for k_ in G[u_][v_] if G[u_][v_][k_]['type']==ETYPE_EXTR]
            assert(len(cand)==1)
            k_ = cand[0]
            if G[u][v][k].get('is_set',True) and G[u_][v_][k_].get('is_set',True):
                to_set.add((u,v,k))
                to_set.add((u_,v_,k_))
            else:
                #bad cycle
                break
            corr_at = {u:v,v:u,v_:u_,u_:v_}
            for x,y in corr_at.items():
                wrong = set([(min(w,x),max(w,x),k) for w in G[x] for k in G[x][w][k] if G[x][w][k]['type']==ETYPE_EXTR and w!=y])
                to_remove.update(wrong)
        if len(to_set.intersection(to_remove))>0:
            #failed cycle
            continue
        violated = [(x,y,z) for x,y,z in to_remove if G[x][y][z].get('is_set',False)]
        if len(violated) > 0:
            continue
        for u,v,k in to_set:
            G[u][v][k]['is_set']=True


def traceback_single_vertex(G, u, trees, x, etype):
    xedges=[]
    aedges=[]
    curr = x
    curr_etype = etype
    while curr != u:
        next = trees[u][curr_etype][curr]
        candidates = [(curr,next,k) for k in G[curr][next] if G[curr][next][k]==curr_etype]
        assert(len(candidates)==1)
        if curr_etype==ETYPE_ADJ:
            aedges.extend(candidates)
        else:
            xedges.extend(candidates)
    return aedges,xedges


def invert_etype(etype):
    return ETYPE_EXTR if etype==ETYPE_ADJ else ETYPE_ADJ


def trace_cycle(G,u,v,x,trees,etype):
    xedges = []
    aedges = []
    ae,xe=traceback_single_vertex(G, u, trees,x, etype)
    aedges.extend(ae)
    xedges.extend(xe)
    ae,xe=traceback_single_vertex(G, v, trees,x, invert_etype(etype))
    return aedges,xedges


def bfs_step(G, trees, frontiers, etype, u,forbidden=None):
    if forbidden is None:
        forbidden=[u]
    newfrontieru = []
    changed = False
    #print(etype,file=stderr)
    for x in frontiers[u]:
        for x_ in G[x]:
            for k in G[x][x_]:
                if G[x][x_][k]['type']==etype and G[x][x_][k].get('is_set',True) and x_ not in trees[u][ETYPE_ADJ] and x_ not in trees[u][ETYPE_EXTR] and x_ not in forbidden:
                    newfrontieru.append(x_)
                    trees[u][etype][x_]=x
                    changed=True
    oldfrontiers = frontiers[u]
    frontiers[u]=newfrontieru
    if changed:
        assert(u in trees[u][ETYPE_EXTR].values())
    return changed


def cycle_greedy_partial_matching(G,sibmap,max_depth=None):
    trees = dict()
    frontiers = dict()
    depths = dict()
    adjacencies = set(((u,v) for u,v,etype in G.edges(data='type') if etype == ETYPE_ADJ and G[u][v].get('is_set',True)))
    for v in G.nodes():
        tree = {ETYPE_ADJ : dict(), ETYPE_EXTR : dict()}
        frontier = [v]
        trees[v]=tree
        frontiers[v]=frontier
        depths[v]=dict()
    changed = True
    etype = ETYPE_ADJ
    depth = 0
    cycles = []
    while changed:
        depth+=1
        etype=invert_etype(etype)
        rm_adj = set()
        changed = False
        for u,v in adjacencies:
            if (u,v) in rm_adj:
                continue
            c=bfs_step(G,trees,frontiers,etype,u,forbidden=[u,v])
            changed=changed or c
            c=bfs_step(G,trees,frontiers,etype,u,forbidden=[u,v])
            changed = changed or c
            for x in frontiers[u]:
                if x in trees[v][invert_etype(etype)]:
                    aedges, xedges = trace_cycle(G,u,v,x,trees,etype)
                    rm_adj.update(aedges)
                    cycles.append(xedges)
        adjacencies.difference_update(rm_adj)
    set_cycles(G,sibmap,cycles)


def try_fix_edges(G,edges,sibmap):
    to_set = set()
    to_remove = set()
    all_edges=set()
    for x,y in edges:
        to_set.add(tuple(sorted((x,y))))
        x_,y_ = sibmap[(x,y)]
        to_set.add(tuple(sorted((x_,y_))))
        #Remove other edges
        rm = [tuple(sorted((w,z))) for w in [x,y,x_,y_] for z in G[w] for k in G[w][z] if G[w][z][k]['type']==ETYPE_EXTR and z not in [x,y,x_,y_]]
        rm_sib = [tuple(sorted(sibmap[e])) for e in rm]
        to_remove.update(rm)
        to_remove.update(rm_sib)
        all_edges.update([tuple(sorted((w,z))) for w in [x,y,x_,y_] for z in G[w] for k in G[w][z] if G[w][z][k]['type']==ETYPE_EXTR])
    try:
        assert(len(to_remove.union(to_set).symmetric_difference(all_edges))==0)
    except AssertionError:
        LOG.debug('To remove: ',to_remove)
        LOG.debug('To set: ',to_set)
        LOG.debug('All: ',all_edges)
        LOG.error('did not consider all edges for all vertices :/')
        raise AssertionError('did not consider all edges for all vertices :/')
    if len(to_set.intersection(to_remove))>0:
        #bad, conflicting edges :/
        return False
    for x,y in to_remove:
        candidates = [k for k in G[x][y] if G[x][y][k]['type']==ETYPE_EXTR]
        assert(len(candidates)==1)
        if  G[x][y][candidates[0]].get('is_set',False):
            return False
    for x,y in to_set:
        candidates = [k for k in G[x][y] if G[x][y][k]['type']==ETYPE_EXTR]
        assert(len(candidates)==1)
        G[x][y][candidates[0]]['is_set']=True
    for x,y in to_remove:
        candidates = [k for k in G[x][y] if G[x][y][k]['type']==ETYPE_EXTR]
        assert(len(candidates)==1)
        G[x][y][candidates[0]]['is_set']=False


def cleanup_tree(G,trees,frontiers,u):
    return
    LOG.debug(f'Change no {u} cleanup {trees[u]}')
    #convert to a parent -> children tree
    pc = {ETYPE_ADJ:dict(),ETYPE_EXTR:dict()}
    if len(trees[u][ETYPE_EXTR])==0:
        #nothing to be done, reset adjacencies and frontiers
        trees[u][ETYPE_ADJ]=dict()
        frontiers[u]=[]
        return
    for etype in [ETYPE_ADJ,ETYPE_EXTR]:
        for child,parent in trees[u][etype].items():
            if not parent in pc[etype]:
                pc[etype][parent]=[]
            pc[etype][parent].append(child)
#    print('ä'*20,file=stderr)
#    print(u,trees[u],file=stderr)
#    print(u,pc,file=stderr)
#    print(u,trees[u][ETYPE_EXTR].values(),file=stderr)
#    print('ä'*20,file=stderr)
    assert(u in trees[u][ETYPE_EXTR].values())
    assert(u in pc[ETYPE_EXTR])
    stack = [(u,ETYPE_EXTR)]
    while len(stack)>0:
        x,etype = stack.pop()
        for child in pc[etype].get(x,[]):
            candidates = [k for k in G[x][child] if G[x][child][k]['type']==etype and G[x][child][k].get('is_set',True)]
            if len(candidates) == 0:
                #dead edge, remove
                remove_dead_branch(trees,frontiers,u,pc,child,invert_etype(etype))
            else:
                stack.append((child,invert_etype(etype)))
#    print('Change no ',u,'cleanup done',trees[u],file=stderr)


def try_fix_cycle(G,sibmap,u,v,trees,frontiers,x,etype):
    edges = []
    aedges = []
    curr = x
    curr_etype=etype
    success = True
    count=0
    while curr != u:
        count+=1
        try:
            x_=trees[u][curr_etype][curr]
        except KeyError:
            success=False
            break
        if curr_etype==ETYPE_EXTR:
            candidates = [k for k in G[curr][x_] if G[curr][x_][k]['type']==ETYPE_EXTR]
            assert(len(candidates)==1)
            if G[curr][x_][candidates[0]].get('is_set',True):
                edges.append((curr,x_))
            else:
                #cycle cannot be used
                cleanup_tree(G,trees,frontiers,u)
                #remove_from_tree(trees,u,x,x_,etype)
                success=False
                break
        else:
            aedges.append((curr,x_))
        curr_etype=invert_etype(curr_etype)
        curr = x_
    curr_etype=invert_etype(etype)
    curr=x
    count=0
    while curr != v:
        count+=1
        try:
            x_=trees[v][curr_etype][curr]
        except KeyError:
            success=False
            break
        if curr_etype==ETYPE_EXTR:
            candidates = [k for k in G[curr][x_] if G[curr][x_][k]['type']==ETYPE_EXTR]
            assert(len(candidates)==1)
            if G[curr][x_][candidates[0]].get('is_set',True):
                edges.append((curr,x_))
            else:
                #cycle cannot be used
                cleanup_tree(G,trees,frontiers,u)
                #remove_from_tree(trees,v,x,x_,invert_etype(etype))
                success=False
                break
        else:
            aedges.append((curr,x_))
        curr_etype=invert_etype(curr_etype)
        curr = x_
    if not success:
        #stop here and do not fix any edges
        return success,[]
    success = try_fix_edges(G,edges,sibmap)
    if not success:
        #TODO: cycle is bad, remove it from tree
        return success,[]
    return success, aedges


def heuristic_partial_matching(G,sibmap,max_iter=None):
    #TODO: Implement max_iter
    trees = dict()
    frontiers = dict()
    #TODO: Filter out useless adjacencies, i.e. those that don't have extremity-edges
    adjacencies = set(((u,v) for u,v,etype in G.edges(data='type') if etype == ETYPE_ADJ and G[u][v].get('is_set',True)))
    for v in G.nodes():
        tree = {ETYPE_ADJ : dict(), ETYPE_EXTR : dict()}
        frontier = [v]
        trees[v]=tree
        frontiers[v]=frontier
    #bfs for each node starting with extremity edge until end of a path is reached or connection to neighboring adjacency found
    changed = True
    etype = ETYPE_ADJ

    while changed:
        etype = invert_etype(etype)
        rm_adj = set()
        changed = False
        for v in G.nodes():
            extren = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==ETYPE_EXTR and G[v][u][k].get('is_set',False)]
            assert(len(extren)<=1)
        for u,v in adjacencies:
            bad=False
            bfs = bfs_step(G, trees, frontiers, etype, u)
            changed = changed or bfs
            bfs = bfs_step(G, trees, frontiers, etype, v)
            changed = changed or bfs
            if len(trees[u][ETYPE_EXTR])==0 or len(trees[v][ETYPE_EXTR]) == 0:
                rm_adj.add((u,v))
                bad=True
            else:
                #TODO: comment out
                for x in (u,v):
                    assert(x in trees[x][ETYPE_EXTR].values())
            #if bf != changed:
            #    print(frontiers[u],file=stderr)
            #    print(trees[v][invert_etype(etype)],file=stderr)
            #    print('-------------',file=stderr)
            #check whether anything special happened at the frontiers, i.e. telomere, indel, cycle
            if not bad:
                for x in frontiers[u]:
                    if x in trees[v][invert_etype(etype)]:
                        #cycle detected
                        for x in (u,v):
                            assert(x in trees[x][ETYPE_EXTR].values())
                        x,aedges=try_fix_cycle(G,sibmap,u,v,trees,frontiers,x,etype)
                        if x:
                            #remove the adjacency edges now dealt with
                            rm_adj.update(aedges)
                            break
                        for x in (u,v):
                            assert(x in trees[x][ETYPE_EXTR].values() or len(trees[x][ETYPE_EXTR])==0)
        adjacencies.difference_update(rm_adj)
    for v in G.nodes():
            extren = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==ETYPE_EXTR and G[v][u][k].get('is_set',False)]
            assert(len(extren)<=1)


def remove_from_tree(trees,u,x,x_,etype):
    curr = x
    curr_etype=etype
    while curr != x_:
        curr=trees[u][curr_etype].pop(curr)
        curr_etype=invert_etype(curr_etype)


def remove_dead_branch(trees,frontiers,u,pc,start,petype):
    dead = {ETYPE_ADJ : [], ETYPE_EXTR : []}
    stack = [(start,petype)]
#    print('*'*20,file=stderr)
#    print('Before',file=stderr)
#    print(start,petype,file=stderr)
#    print(trees[u],file=stderr)
#    print(pc,file=stderr)
#    print('*'*20,file=stderr)
    while len(stack)>0:
        x,etype = stack.pop()
        #print('Stacklen', len(stack),file=stderr)
        dead[invert_etype(etype)].append(x)
        for child in pc[etype].get(x,[]):
            stack.append((child,invert_etype(etype)))
    for etype, deads in dead.items():
        frontiers_new = set(frontiers[u]).difference(set(deads))
        frontiers[u]=frontiers_new
        for d in deads:
#            print(u,etype,d,file=stderr)
            trees[u][etype].pop(d)
#    print('*'*20,file=stderr)
#    print('After',file=stderr)
#    print(start,petype,file=stderr)
#    print(trees[u],file=stderr)
#    print(pc,file=stderr)
#    print('*'*20,file=stderr)


def fits_ptype(G,genomes,mx,r):
    if r.isupper() != (G.nodes[mx]['type']==VTYPE_CAP):
        return False
    if (r.upper()=='A') != (get_genome(G,mx) == genomes[0]):
        return False
    return True


def sol_from_decomposition(graphs,circ_singletons,alpha,out):
    #TODO: Implement
    var_map = dict()
    wsum = 0.0
    fsum = 0
    for tree_edge, ((child, parent), G) in enumerate(sorted(graphs.items())):
        genomes = [child,parent]
        local = G.copy()
        counts = {'2n':0,'w':0.0,'c':0,'ab':0,'AB':0,'Aa':0,'Ab':0,'Ba':0,'Bb':0,'s':0}
        for u,v,k,d in G.edges(keys=True,data=True):
            xvar = 'x{sep}{te}{sep}{e}'.format(sep=SEP,te=tree_edge,e=d['id'])
            if not d.get('is_set',False) or d['type']==ETYPE_ID:
                local.remove_edge(u,v,key=k)
            var_map[xvar]= 1 if d.get('is_set',False) else 0
            if d['type']==ETYPE_ADJ:
                ancvar = 'a{sep}{e}'.format(sep=SEP,e=G[u][v][k]['anc'])
                if not ancvar in var_map:
                    var_map[ancvar] =  1 if d.get('is_set',False) else 0
                else:
                    assert(var_map[ancvar] ==  (1 if d.get('is_set',False) else 0))
                if d.get('is_set',False):
                    counts['w']+=d['weight']
            if d['type']==ETYPE_EXTR and d.get('is_set',False):
                counts['2n']+=1
            assert(local.has_edge(u,v,key=k)==(var_map[xvar]==1) or d['type']==ETYPE_ID)
        if (child,parent) in circ_singletons:
            for j,cs in enumerate(circ_singletons[(child, parent)].values()):
                is_set = 1
                for data in cs:
                    if var_map['x{sep}{te}{sep}{e}'.format(sep=SEP,te=tree_edge,e=data['id'])]==0:
                        is_set=0
                        break
                csvar='rms{sep}{te}{sep}{j}'.format(sep=SEP,te=tree_edge,j=j)
                var_map[csvar]=is_set
                counts['s']+=is_set
        else:
            #TODO: implement other circ sing counting
            for x, data in G.nodes(data=True):
                if data.get('cscandidate',False):
                    counts['s']+=1
            #local=G.copy()
            #for u,v,k,data in G.edges(keys=True,data=True):
            #    if data['type']==ETYPE_EXTR or not data.get('is_set',False):
            #        local.remove_edge(u,v,k)
            #for comp in nx.connected_components(local):
            #    p_ends = []
            #    for v in comp:
            #        dg = local.degree[v]
            #        assert(dg <=2)
            #        if dg==1:
            #            p_ends.append(v)
            #    assert(len(p_ends)<=2)
            #    if len(p_ends)==1:
            #        continue
            #    if len(p_ends)==0:
            #        for x in nx.edge_dfs(list(comp)[0]):
            #            pass


        for u,data in G.nodes(data=True):
            ancvar = 'g{sep}{anc}'.format(sep=SEP,anc=data['anc'])
            acands = [v for v in G[u] for k in G[u][v] if G[u][v][k]['type']==ETYPE_ADJ and G[u][v][k].get('is_set',False)]
            ecands = [v for v in G[u] for k in G[u][v] if G[u][v][k]['type']==ETYPE_EXTR and G[u][v][k].get('is_set',False)]
            icands = [v for v in G[u] for k in G[u][v] if G[u][v][k]['type']==ETYPE_ID and G[u][v][k].get('is_set',False)]
            assert(len(acands)>=len(ecands)+len(icands))
            assert(len(acands)<=1)
            if len(acands)==0:
                ugnm = get_genome(G,u)
                var_map['z{sep}{te}{sep}{v}'.format(sep=SEP,te=tree_edge,v=u)]=0
                var_map['y{sep}{te}{sep}{v}'.format(sep=SEP,te=tree_edge,v=u)]=1
                var_map['l{sep}{te}{sep}{v}'.format(sep=SEP,te=tree_edge,v=u)]=0 if ugnm==genomes[0] else 1
                if data['type']==VTYPE_CAP:
                    if ugnm==genomes[0]:
                        reps = ['Ab','Aa','AB']
                    else:
                        reps = ['Bb','Ba']
                    for rep in reps:
                        var_map['r{r}{sep}{te}{sep}{v}'.format(r=rep,sep=SEP,te=tree_edge,v=u)]=0
                elif ugnm ==genomes[0]:
                    var_map['rab{sep}{te}{sep}{v}'.format(r=rep,sep=SEP,te=tree_edge,v=u)]=0
                    var_map['rc{sep}{te}{sep}{v}'.format(r=rep,sep=SEP,te=tree_edge,v=u)]=0
            if not ancvar in var_map:
                var_map[ancvar]=len(acands)
            else:
                assert(var_map[ancvar]==len(acands))
        for comp in nx.connected_components(local):
            path_ends=[]
            min_id = None
            for v in comp:
                assert(local.degree(v)<=2)
                if local.degree(v)<2:
                    path_ends.append(v)
                min_id = v if min_id is None else min(min_id,v)
            #TODO: single telomeres
            if len(path_ends)==1:
                v=list(comp)[0]
                assert(G.nodes[path_ends[0]]['type']==VTYPE_CAP or G.nodes[path_ends[0]]['is_set']==False)
                continue
            assert(len(path_ends)==2 or len(path_ends)==0)
            is_positive=True
            for v in path_ends:
                if G.nodes[v]['type']!=VTYPE_CAP:
                    is_positive=False
            if is_positive:
                for v in comp:
                    var_map['y{sep}{te}{sep}{v}'.format(v=v,sep=SEP,te=tree_edge)]=min_id
                    var_map['z{sep}{te}{sep}{v}'.format(v=v,sep=SEP,te=tree_edge)]= 0 if v!=min_id else 1
            else:
                for v in comp:
                    var_map['y{sep}{te}{sep}{v}'.format(v=v,sep=SEP,te=tree_edge)]=0
                    var_map['z{sep}{te}{sep}{v}'.format(v=v,sep=SEP,te=tree_edge)]=0
            if len(path_ends)==0:
                #cycle
                counts['c']+=1
                for v in comp:
                    var_map['l{sep}{te}{sep}{v}'.format(v=v,sep=SEP,te=tree_edge)]=0
                    if get_genome(G,v) == genomes[0]:
                        var_map['rc{sep}{te}{sep}{v}'.format(sep=SEP,te=tree_edge,v=v)] = 0 if v!=min_id else 1
                        var_map['rab{sep}{te}{sep}{v}'.format(sep=SEP,v=v,te=tree_edge)]=0
            else:
                mn = min(path_ends)
                mx = max(path_ends)
                for v in comp:
                    if get_genome(G,v)==genomes[0] and G.nodes[v]['type']!=VTYPE_CAP:
                        #not a cycle
                        var_map['rc{sep}{te}{sep}{v}'.format(sep=SEP,te=tree_edge,v=v)]=0
                if get_genome(G,mn) == get_genome(G,mx):
                #same end, so can just pick this as the genome
                    for v in comp:
                        var_map['l{sep}{te}{sep}{v}'.format(v=v,sep=SEP,te=tree_edge)]= 0 if get_genome(G,mx) == genomes[0] else 1
                else:
                    mxg = 0 if get_genome(G,mx)==genomes[0] else 1
                    mng = 0 if get_genome(G,mn)==genomes[0] else 1
                    for v in comp:
                        var_map['l{sep}{te}{sep}{v}'.format(v=v,sep=SEP,te=tree_edge)]= mxg if v!=mn else mng
                for v in comp:
                    if v in [mx,mn] or get_genome(G,v)!=genomes[0]:
                        continue
                    var_map['rab{sep}{te}{sep}{v}'.format(sep=SEP,v=v,te=tree_edge)]=0
                #set reporting variables at ends
                #TODO: What about even paths?
                if G.nodes[mn]['type']!=VTYPE_CAP:
                    #must be type ab, because caps generally have lower ids
                    assert(G.nodes[mx]['type']!=VTYPE_CAP)
                    assert(get_genome(G,mn)==genomes[0] or get_genome(G,mn)==get_genome(G,mx))
                    if get_genome(G,mn)==genomes[0]:
                        if get_genome(G,mn)!=get_genome(G,mx):
                            var_map['rab{sep}{te}{sep}{v}'.format(sep=SEP,v=mn,te=tree_edge)]=1
                            counts['ab']+=1
                        else:
                            var_map['rab{sep}{te}{sep}{v}'.format(sep=SEP,v=mn,te=tree_edge)]=0
                            var_map['rab{sep}{te}{sep}{v}'.format(sep=SEP,v=mx,te=tree_edge)]=0
                    #do not uncomment this var_map['rab{sep}{te}{sep}{v}'.format(sep=SEP,v=mx,te=tree_edge)]=0
                else:
                    assert(min_id==mn)
                    if G.nodes[mx]['type']!=VTYPE_CAP:
                        #only one var to set
                        if get_genome(G,mx)==genomes[0]:
                            var_map['rab{sep}{te}{sep}{v}'.format(sep=SEP,v=mx,te=tree_edge)]=0
                    else:
                        if get_genome(G,mx)==genomes[0]:
                            reps = ['Ab','Aa','AB']
                        else:
                            reps = ['Ba','Bb']
                        for r in reps:
                            var_map['r{r}{sep}{te}{sep}{v}'.format(r=r,sep=SEP,v=mx,te=tree_edge)]=0
                    if get_genome(G,mn)==genomes[0]:
                        reps = ['Ab','Aa','AB']
                    else:
                        reps = ['Ba','Bb']
                    for r in reps:
                        var_map['r{r}{sep}{te}{sep}{v}'.format(r=r,sep=SEP,v=mn,te=tree_edge)]=1 if fits_ptype(G,genomes,mx,r[1]) else 0
                        counts[r]+=var_map['r{r}{sep}{te}{sep}{v}'.format(r=r,sep=SEP,v=mn,te=tree_edge)]
            if len(path_ends)>0:
                mnreportv = ['r{tp}{sep}{te}{sep}{v}'.format(sep=SEP,te=tree_edge,v=mn,tp=tp) for tp in ['AB','Aa','Ab','Ba','Bb','ab']]
                mnreportvl = [var_map.get(rv,None) for rv in mnreportv]
                mnl = var_map['l{sep}{te}{sep}{v}'.format(sep=SEP,te=tree_edge,v=mn)]
                mxl = var_map['l{sep}{te}{sep}{v}'.format(sep=SEP,te=tree_edge,v=mx)]
                assert(1 in mnreportvl or (mxl == mnl and (G.nodes[mn]['type'],get_genome(G,mn)) == (G.nodes[mx]['type'],get_genome(G,mx))))
        n = int(ceil(counts['2n']/2))
        var_map['n{sep}{te}'.format(sep=SEP,te=tree_edge)]=n
        var_map['c{sep}{te}'.format(sep=SEP,te=tree_edge)]=counts['c']
        var_map['w{sep}{te}'.format(sep=SEP,te=tree_edge)]=counts['w']
        var_map['s{sep}{te}'.format(sep=SEP,te=tree_edge)]=counts['s']
        wsum+=counts['w']
        for r in ['ab','AB','Aa','Ab','Ba','Bb']:
            var_map['p{r}{sep}{te}'.format(r=r,sep=SEP,te=tree_edge)]=counts[r]
        pABa = max(counts['Aa'],counts['Ba'])
        pABb = max(counts['Ab'],counts['Bb'])
        var_map['pABa{sep}{te}'.format(sep=SEP,te=tree_edge)]=pABa
        var_map['pABb{sep}{te}'.format(sep=SEP,te=tree_edge)]=pABb
        q = int(ceil((counts['ab']+pABa+pABb - counts['AB'])/2))
        var_map['q{sep}{te}'.format(sep=SEP,te=tree_edge)]=q
        #TODO: Circular singletons!
        f = n - counts['c'] + q + counts['s']
        var_map['f{sep}{te}'.format(sep=SEP,te=tree_edge)]=f
        fsum+=f
    obj = alpha*fsum + (alpha-1)*wsum
    for tree_edge, ((child, parent), G) in enumerate(sorted(graphs.items())):
        genomes = [child,parent]
        for v,data in G.nodes(data=True):
            if get_genome(G,v)!=genomes[0] or data['type']==VTYPE_CAP:
                continue
            assert('rab{sep}{te}{sep}{v}'.format(sep=SEP,te=tree_edge,v=v) in var_map)
    print('# Objective value = {}'.format(obj),file=out)
    for vr,vl in var_map.items():
        print(vr,vl,file=out)


def greedy_extend_max_match(local,mwmtch):
    mwmap = dict()
    for x,y in mwmtch:
        mwmap[x]=y
        mwmap[y]=x
    for v,data in local.nodes(data=True):
        if data['type']==VTYPE_CAP:
            #caps don't need to be matched
            continue
        if v in mwmap:
            continue
        #print('Not in Map',file=stderr)
        candidates = [(local[v][u]['weight'],u) for u in local[v] if u not in mwmap]
        assert(len(candidates)>0)
        u = max(candidates,key=lambda x:x[0])[1]
        mwmap[v]=u
        mwmap[u]=v
    #print(mwmap,file=stderr)
    return list(set((tuple(sorted(x)) for x in mwmap.items())))


def extract_local_graph(genome, G):
    local = G.copy()
    for v in G.nodes():
        if get_genome(G,v) != genome:
            local.remove_node(v)
    return local


def get_max_match(genome, G):
    medges = []
    local = extract_local_graph(genome, G)
    for u,v,tp in list(local.edges(data='type')):
        if tp!=ETYPE_ADJ:
            local.remove_edge(u,v)
    for v in list(G.nodes()):
        if not G.nodes[v].get('is_set',True):
            G.remove_node(v)
    telomeres = [u for u,data in G.nodes(data=True) if data['type']==VTYPE_CAP]
    telomeres = set(telomeres)
    #set every edge next to a regular vertex to more than weight 1
    weights = [w for u,v,w in local.edges(data='weight')]
    max_weight = max(weights)
    min_weight = min(weights)
    number_es = len(local.nodes())
    for u,v,k in local.edges(keys=True):
        w=(local[u][v][k]['weight']-min_weight)/(max_weight+1)
        if u in telomeres or v in telomeres:
            mw = number_es+w
        else:
            mw = 2*number_es+w

        local[u][v][k]['mod_weight']=mw
    #this works since adjacency edges generally do not double
    local = nx.Graph(local)
    #print('\n'.join([str(k) for k in sorted(G.nodes(data=True))]),file=stderr)
    LOG.info('calculating maximum weight matching')
    mwmtch = list(nx.max_weight_matching(local,weight='mod_weight'))
    LOG.info('maximum weight matching done')
    #print(mwmtch,file=stderr)
    #check matching
    matched = set([v for v,_ in mwmtch]+[v for _,v in mwmtch])
    for v,tp in local.nodes(data='type'):
        if tp == VTYPE_CAP:
            continue
        assert(v in matched)
    #mwmtch=greedy_extend_max_match(local,mwmtch)
    #print(mwmtch,file=stderr)
    return [(G.nodes[x]['anc'],G.nodes[y]['anc']) for x,y in mwmtch]


def get_anc_map(g):
    return dict(((anc,v) for v,anc in g.nodes(data='anc')))


def set_based_on_matching(G,mtch,genome):
    local = extract_local_graph(genome, G)
    ancmap = get_anc_map(local)
    for u_,v_ in mtch:
        u,v = ancmap[u_],ancmap[v_]
        candidates = [k for k in G[u][v] if G[u][v][k]['type']==ETYPE_ADJ]
        assert(len(candidates)==1)
        G[u][v][candidates[0]]['is_set']=True
    for u,v,k,data in local.edges(keys=True,data=True):
        if data['type']!=ETYPE_ADJ:
            #skip indel edges
            continue
        if not 'is_set' in G[u][v][k]:
            G[u][v][k]['is_set']=False


def annotate_superfluous(G,sup_ancs):
    anc_to_v = dict()
    for v,anc in G.nodes(data='anc'):
        anc_to_v[anc]=v
    for anc in sup_ancs:
        v=anc_to_v[anc]
        G.nodes[v]['is_set']=False


def find_other(G,sibmap,v):
    idec = [u for u in G[v] for k in G[v][u] if G[v][u][k]['type']==ETYPE_ID]
    assert(len(idec)<=1)
    if len(idec)==1:
        return idec[0]
    extc = [u for u in G[v] for k in G[v][u] if G[v][u][k]['type']==ETYPE_EXTR]
    assert(len(extc)>=1)
    u=extc[0]
    x,y=sibmap[(v,u)]
    if get_genome(G,x)==get_genome(G,v):
        return x
    else:
        return y


def get_sup_families(G,not_full_fams,gnm):
    nff_markers = dict()
    #if the matching does not allow for a good node
    for v,data in G.nodes(data=True):
        f = nodeGetFamily(G,v)
        if data['type']==VTYPE_CAP or get_genome(G,v)!=gnm:
            continue
        if data['id'][1][1]!=EXTR_TAIL:
            continue
        if f in not_full_fams[gnm]:
            if f not in nff_markers:
                nff_markers[f]=set()
            nff_markers[f].add(v)
    return nff_markers


def get_superfluous_ancs(G,sibmap,not_full_fams,gnm):
    sups = get_sup_families(G,not_full_fams,gnm)
    max_weights = dict()
    ancs = set()
    for f,heads in sups.items():
        for h in heads:
            t = find_other(G,sibmap,h)
            mt = max([G[t][u][k]['weight'] for u in G[t] for k in G[t][u] if G[t][u][k]['type']==ETYPE_ADJ])
            mh = max([G[h][u][k]['weight'] for u in G[h] for k in G[h][u] if G[h][u][k]['type']==ETYPE_ADJ])
            max_weights[(t,h)]=mt+mh
        worst=sorted(list(max_weights.keys()),key=lambda x: max_weights[x])
        to_remove=not_full_fams[gnm][f]
        for t,h in worst[:to_remove]:
            ancs.add(G.nodes[t]['anc'])
            ancs.add(G.nodes[h]['anc'])
    return ancs


def create_sibmap(G,siblings):
    id_to_edges = dict()
    sibmap=dict()
    for u,v,data in G.edges(data=True):
        if data['type']==ETYPE_EXTR:
            id_to_edges[data['id']]=(u,v)
    #print(siblings,file=stderr)
    for siba,sibb in siblings:
        if not siba in id_to_edges or not sibb in id_to_edges:
            continue
        x,y=id_to_edges[siba]
        w,z =id_to_edges[sibb]
        for u,v in [(x,y),(y,x)]:
            sibmap[(u,v)]=(w,z)
        for u,v in [(w,z),(z,w)]:
            sibmap[(u,v)]=(x,y)
    return sibmap


def create_max_weight(graphs,not_full_fams,siblings):
    root_candidates = dict()
    not_root = set()
    matchings = dict()
    superfluous_anc = dict()
    for tree_edge, ((child, parent), G) in enumerate(sorted(graphs.items())):
        sibmap=create_sibmap(G,siblings[(child,parent)])
        root_candidates[parent]=child
        not_root.add(child)
        superfluous_anc[child] = get_superfluous_ancs(G,sibmap,not_full_fams,child)
        annotate_superfluous(G,superfluous_anc[child])
        matchings[child]=get_max_match(child, G)
    for x in not_root:
        root_candidates.pop(x,None)
    assert(len(root_candidates)==1)
    parent,child = list(root_candidates.items())[0]
    G=graphs[(child,parent)]
    sibmap=create_sibmap(G,siblings[(child,parent)])
    superfluous_anc[parent] = get_superfluous_ancs(G,sibmap,not_full_fams,parent)
    annotate_superfluous(G,superfluous_anc[parent])
    matchings[parent]=get_max_match(parent,G)
    for tree_edge, ((child, parent), G) in enumerate(sorted(graphs.items())):
        annotate_superfluous(G,superfluous_anc[child])
        set_based_on_matching(G,matchings[child],child)
        annotate_superfluous(G,superfluous_anc[parent])
        set_based_on_matching(G,matchings[parent],parent)


def warm_start_decomposition(graphs,fam_bounds,families,siblings):
    not_full_fams = getNotFullFamilies(fam_bounds,families)
    create_max_weight(graphs,not_full_fams,siblings)
    for tree_edge, ((child, parent), G) in enumerate(sorted(graphs.items())):
        sibmap=create_sibmap(G,siblings[(child,parent)])
        #cycle_greedy_partial_matching(G,sibmap)
        heuristic_partial_matching(G,sibmap)
        greedy_matching_completion(G,sibmap)


def add_weight_tels(G,add_telweight):
    for u,v,k,data in G.edges(keys=True,data=True):
        if data['type']==ETYPE_ADJ and VTYPE_CAP in [G.nodes[x]['type'] for x in [u,v]]:
            G[u][v][k]['weight']=G[u][v][k]['weight']+add_telweight

