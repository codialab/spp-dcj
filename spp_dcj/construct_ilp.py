# import from built-in packages
from argparse import FileType
from sys import stdout
from itertools import chain
from functools import reduce
import logging

# import from third-party packages
import networkx as nx

# import from own packages
import spp_dcj.data_utils as du


#
# global variables
#

ADJ_TRUST_THRESHOLD = 0.9

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

#
# ILP OBJECTIVE
#

def objective(graphs, alpha, out):
    out.write('minimize ')
    terms = ["{ialph} w{sep}{te} + {alph} f{sep}{te}".format(
        alph=alpha,ialph=alpha-1,sep=du.SEP,te=tree_edge)
        for tree_edge, _ in enumerate(sorted(graphs.items()))]
    print(" + ".join(terms),file=out)



def child_parent_tree(graphs):
    cp_tree = dict()
    for i, ((child, parent), G) in enumerate(sorted(graphs.items())):
        assert(child not in cp_tree)
        cp_tree[child]=(parent,i)
    return cp_tree


def trace_to_root(cp_tree,x):
    trace = dict()
    curr=x
    while curr in cp_tree:
        next, tree_edge = cp_tree[curr]
        trace[curr]=(next,tree_edge)
        curr = next
    trace[curr]=(None,None)
    return trace

def lca_treeedges_cp_tree(cp_tree,a,b):
    #trace a to root
    trace_a = trace_to_root(cp_tree,a)
    curr = b
    edges = []
    while curr not in trace_a:
        next,tree_edge = cp_tree[curr]
        edges.append(tree_edge)
        curr = next
    end = curr
    curr = a
    while curr!=end:
        next,tree_edge = trace_a[curr]
        edges.append(tree_edge)
        curr=next
    return edges

#
# ILP CONSTRAINTS
#

def constraints(graphs, siblings,circ_singletons, out,lower_bound_mat={}):
    out.write('subject to\n')
    for i, ((child, parent), G) in enumerate(sorted(graphs.items())):

        LOG.info(('writing constraints for relational diagram of {} and ' + \
                '{}').format(child, parent))
        genomes = [child,parent]
        global_constraints(G,out)
        loose_constraints(G,i,genomes,out)
        slm_constraints(G,i,siblings,out)
        regular_reporting(G,i,genomes,out)
        pcap_reporting(G,i,genomes,out)
        if (child,parent) in circ_singletons:
            c18(G,i,out,candidates=circ_singletons[(child,parent)])
            manual_cs_constraints(circ_singletons[(child,parent)],i,out)
        else:
            c18(G,i,out)
            cs_constraints(G,i,out,max_circ_len=len(G.nodes())/2)
    LOG.info("Writing lower bound constraints.")
    cp_tree=child_parent_tree(graphs)
    print(cp_tree,file=stderr)
    for a in lower_bound_mat:
        for b in lower_bound_mat[a]:
            lb = lower_bound_mat[a][b]
            tree_edges=lca_treeedges_cp_tree(cp_tree,a,b)
            if len(tree_edges)==0:
                continue
            fvars = ['f{sep}{te}'.format(sep=du.SEP,te=te) for te in tree_edges]
            fsum = ' + '.join(fvars)
            print("{fsum} >= {lb}".format(fsum=fsum,lb=lb))


    out.write('\n')

def get_gene_extremities(G):
    ids_to_v = {data['id'] : v for v,data in G.nodes(data=True)}
    l = []
    for v,data in G.nodes(data=True):
        if data['type'] == du.VTYPE_CAP:
            continue
        l.append((v,ids_to_v[du.complement_id(data['id'])]))
    return l


#TODO: Integrate multiplicities!
def c01(G,out):
    for v,data in G.nodes(data=True):
        if data['type']==du.VTYPE_CAP:
            continue
        print("g{sep}{v} = 1".format(sep=du.SEP,v=data['anc']),file=out)

#Old,removed
#def c02(G,out):
#    for u,v in get_gene_extremities(G):
#        print("g{sep}{u} - g{sep}{v} = 0".format(u=G.nodes[u]['anc'],v=G.nodes[v]['anc'],sep=du.SEP),file=out)

def c02(G,out):
    for v,data in G.nodes(data=True):
        aes = [G[v][u][i]['anc'] for u in G[v] for i in G[v][u] if G[v][u][i]['type']==du.ETYPE_ADJ]
        sm = ' + '.join(['a{sep}{e}'.format(sep=du.SEP,e=e) for e in aes])
        print("{sm} - g{sep}{v} = 0".format(sm=sm,sep=du.SEP,v=G.nodes[v]['anc']),file=out)


def global_constraints(G,out):
    for c in [c01,c02]:
        c(G,out)

def c03(G,tree_edge,out):
    ws = [ '{w} x{sep}{te}{sep}{e}'.format(w=data['weight'],te=tree_edge,e=data['id'],sep=du.SEP) for u,v,data in G.edges(data=True) if data['type']==du.ETYPE_ADJ]
    print("{s} - w{sep}{treeedge} = 0".format(treeedge=tree_edge,s=' + '.join(ws),sep=du.SEP),file=out)


def c04(G,tree_edge,out):
    print("n{sep}{te} - c{sep}{te} + q{sep}{te} + s{sep}{te} - f{sep}{te} = 0".format(sep=du.SEP,te=tree_edge),file=out)

def c05(G,tree_edge,out):
    xs = ['0.5 x{sep}{te}{sep}{e}'.format(te=tree_edge,sep=du.SEP,e=data['id']) for u,v,data in G.edges(data=True) if data['type']==du.ETYPE_EXTR]
    print('{s} - n{sep}{te} = 0'.format(s=' + '.join(xs),sep=du.SEP,te=tree_edge),file=out)

def c06(G,tree_edge,genomes,out):
    cs = ['rc{sep}{te}{sep}{v}'.format(te=tree_edge,sep=du.SEP,v=v) for v,data in G.nodes(data=True) if data['type'] != du.VTYPE_CAP and get_genome(G,v)==genomes[0]]
    print('{s} - c{sep}{te} = 0'.format(s=' + '.join(cs),sep=du.SEP,te=tree_edge),file=out)

def c07(G,tree_edge,out):
    print("pab{sep}{te} + pABa{sep}{te} + pABb{sep}{te} - pAB{sep}{te} - 2 q{sep}{te} <= 0".format(sep=du.SEP,te=tree_edge),file=out)


def summary_constraint(pathtype,predicate,G,tree_edge,out):
    sm = ['r{ptype}{sep}{te}{sep}{v}'.format(te=tree_edge,ptype=pathtype,sep=du.SEP,v=v) for v in G.nodes() if predicate(G,v)]
    print("{s} - p{ptype}{sep}{te} = 0".format(s=' + '.join(sm),ptype=pathtype,sep=du.SEP,te=tree_edge),file=out)

c08 = lambda G,tree_edge,genomes,out: summary_constraint('ab',lambda G,v: get_genome(G,v)==genomes[0] and G.nodes[v]['type']==du.VTYPE_EXTR,G,tree_edge,out)
c09 = lambda G,tree_edge,genomes,out: summary_constraint('Ab',lambda G,v: get_genome(G,v)==genomes[0] and G.nodes[v]['type']==du.VTYPE_CAP,G,tree_edge,out)    
c10 = lambda G,tree_edge,genomes,out: summary_constraint('Ba',lambda G,v: get_genome(G,v)==genomes[1] and G.nodes[v]['type']==du.VTYPE_CAP,G,tree_edge,out)    
c11 = lambda G,tree_edge,genomes,out: summary_constraint('Aa',lambda G,v: get_genome(G,v)==genomes[0] and G.nodes[v]['type']==du.VTYPE_CAP,G,tree_edge,out)    
c12 = lambda G,tree_edge,genomes,out: summary_constraint('Bb',lambda G,v: get_genome(G,v)==genomes[1] and G.nodes[v]['type']==du.VTYPE_CAP,G,tree_edge,out)    

def c13toc16(G,tree_edge,out):
    for indel in "ab":
        for tel in "AB":
            print("pAB{ind}{sep}{te} - p{tel}{ind}{sep}{te} >= 0".format(ind=indel,tel=tel,te=tree_edge,sep=du.SEP),file=out)
c17 = lambda G,tree_edge,genomes,out: summary_constraint('AB',lambda G,v: get_genome(G,v)==genomes[0] and G.nodes[v]['type']==du.VTYPE_CAP,G,tree_edge,out)

def c18(G,tree_edge,out,candidates=None):
    if candidates is None:
        sm = ['rs{sep}{te}{sep}{v}'.format(te=tree_edge,sep=du.SEP,v=v) for v,data in G.nodes(data=True) if data['type'] != du.VTYPE_CAP and data.get('cscandidate',False)]
        
    else:
        LOG.info("manual candidates: {}".format(candidates))
        sm = ['rms{sep}{te}{sep}{j}'.format(sep=du.SEP,te=tree_edge,j=j) for j in range(len(candidates))]
        LOG.info("sm: {}".format(sm))
    print("{s} - s{sep}{te} = 0".format(s=' + '.join(sm),sep=du.SEP,te=tree_edge),file=out)

def c19(G,tree_edge,out):
    for v,data in G.nodes(data=True):
        if data['type']==du.VTYPE_CAP:
            continue
        ees = [G[v][u][i]['id'] for u in G[v] for i in G[v][u] if G[v][u][i]['type'] in [du.ETYPE_EXTR,du.ETYPE_ID]]
        sm = ' + '.join(['x{sep}{te}{sep}{e}'.format(sep=du.SEP,te=tree_edge,e=e) for e in ees])
        print("{sm} - g{sep}{anc} = 0".format(sm=sm,sep=du.SEP,anc = G.nodes[v]['anc']),file=out)

def c20(G,tree_edge,out):
    for u,v,data in G.edges(data=True):
        if data['type'] == du.ETYPE_ADJ:
            print("a{sep}{anc} - x{sep}{te}{sep}{e} = 0".format(sep=du.SEP,te=tree_edge,e=data['id'],anc=data['anc']),file=out)

def c21(G,tree_edge,out):
    for v,data in G.nodes(data=True):
        print("z{sep}{te}{sep}{v} - g{sep}{anc} <= 0".format(sep=du.SEP,te=tree_edge,v=v,anc=data['anc']),file=out)


def loose_constraints(G,tree_edge,genomes,out):
    for c in [c03,c04,c05,c07,c13toc16,c19,c20,c21]:
        c(G,tree_edge,out)
    for c in [c06,c08,c09,c10,c11,c12,c17]:
        c(G,tree_edge,genomes,out)

def cslm25(sibs,tree_edge,out):
    for eid,did in sibs:
        print("x{sep}{te}{sep}{eid} - x{sep}{te}{sep}{did} = 0".format(sep=du.SEP,eid=eid,did=did,te=tree_edge),file=out)

def cslm26(G,tree_edge,out):
    for u_,v_,data in G.edges(data=True):
        eid = data['id']
        for u,v in [(u_,v_),(v_,u_)]:
            if data['type']==du.ETYPE_ID:
                print("y{sep}{te}{sep}{u} + {u} x{sep}{te}{sep}{eid} <= {u}".format(sep=du.SEP,eid=eid,u=u,te=tree_edge),file=out)
            else:
                print("y{sep}{te}{sep}{v} - {u} x{sep}{te}{sep}{eid} - y{sep}{te}{sep}{u} >= -{u}".format(sep=du.SEP,eid=eid,v=v,u=u,te=tree_edge),file=out)

def cslm27(G,tree_edge,out):
    for v in G.nodes():
        print("{v} z{sep}{te}{sep}{v} - y{sep}{te}{sep}{v} <= 0".format(sep=du.SEP,v=v,te=tree_edge),file=out)

def slm_constraints(G,tree_edge,sibs,out):
    """
    Shao-Lin-Moret Constraints as defined in Table 1.
    """
    cslm25(sibs,tree_edge,out)
    cslm26(G,tree_edge,out)
    cslm27(G,tree_edge,out)






def crv28(G,tree_edge,genomes,out):
    for x,y, data in G.edges(data=True):
        if data['type'] == du.ETYPE_ID:
            for v in [x,y]:
                if get_genome(G, v) == genomes[0]:
                    print("l{sep}{te}{sep}{v} +  x{sep}{te}{sep}{e} <= 1".format(sep=du.SEP,te=tree_edge,v=v,e=data['id']),file=out)
                else:
                    print("l{sep}{te}{sep}{v} -  x{sep}{te}{sep}{e} >= 0".format(sep=du.SEP,te=tree_edge,v=v,e=data['id']),file=out)

def get_genome(G, v):
    return G.nodes[v]['id'][0]


def crv29(G,tree_edge,genomes,out):
    for x,y, data in G.edges(data=True):
        if data['type'] == du.ETYPE_ID or du.VTYPE_CAP in [G.nodes[v]['type'] for v in [x,y]]:
            continue

        for u,v in [(x,y),(y,x)]:
            if get_genome(G,v) == genomes[0]  and data['type'] == du.ETYPE_ADJ:
                print("l{sep}{te}{sep}{v} - l{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} - rab{sep}{te}{sep}{v} - rab{sep}{te}{sep}{u} <= 1".format(
                    sep=du.SEP,
                    te=tree_edge,
                    e=data['id'],
                    u=u,
                    v=v
                ), file=out)
            else:
                print("l{sep}{te}{sep}{v} - l{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} <= 1".format(
                    sep=du.SEP,
                    te=tree_edge,
                    e=data['id'],
                    u=u,
                    v=v
                ), file=out)

def crv30(G,tree_edge,genomes,out):
    for v,data in G.nodes(data=True):
        if data['type']==du.VTYPE_CAP or get_genome(G,v)!=genomes[0]:
            continue
        print('rc{sep}{te}{sep}{v} - z{sep}{te}{sep}{v} <= 0'.format(
            sep=du.SEP,
            te=tree_edge,
            v=v
        ), file=out)

def crv31(G,tree_edge,genomes,out):
    for x, y, data in G.edges(data=True):
        if data['type'] != du.ETYPE_ID or get_genome(G,x) != genomes[0]:
            continue
        for v in [x,y]:
            print("rab{sep}{te}{sep}{v} - x{sep}{te}{sep}{e} <= 0".format(
                sep=du.SEP,
                te=tree_edge,
                v=v,
                e=data['id']
            ),file=out)

def regular_reporting(G,tree_edge,genomes,out):
    for c in [crv28,crv29,crv30,crv31]:
        c(G,tree_edge,genomes,out)


def crt32(G,tree_edge,genomes,out):
    for v,data in G.nodes(data=True):
        if data['type']!= du.VTYPE_CAP:
            continue
        print("l{sep}{te}{sep}{v} = {gid}".format(
            sep=du.SEP,
            te=tree_edge,
            v=v,
            gid=0 if get_genome(G,v)==genomes[0] else 1
        ),file=out)

def crt33(G,tree_edge,genomes,out):
    for x,y,data in G.edges(data=True):
        if data['type']!=du.ETYPE_ADJ or not du.VTYPE_CAP in [G.nodes[v]['type'] for v in [x,y]]:
            continue
        if get_genome(G,x) == genomes[0]:
            if G.nodes[x]['type'] == du.VTYPE_CAP:
                v,u = x,y
            else:
                v,u = y,x
            print("l{sep}{te}{sep}{u} -  l{sep}{te}{sep}{v} - rAB{sep}{te}{sep}{v} - rAb{sep}{te}{sep}{v} + x{sep}{te}{sep}{e} <= 1".format(
                sep=du.SEP,
                te=tree_edge,
                u=u,
                v=v,
                e=data['id']
            ),file=out)
        else:
            if G.nodes[x]['type'] == du.VTYPE_CAP:
                u,v = x,y
            else:
                u,v = y,x
            print("l{sep}{te}{sep}{u} -  l{sep}{te}{sep}{v} - rBa{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} <= 1".format(
                sep=du.SEP,
                te=tree_edge,
                u=u,
                v=v,
                e=data['id']
            ),file=out)

def crt34(G,tree_edge,genomes,out):
    for v,data in G.nodes(data=True):
        if get_genome(G,v)==genomes[0] and data['type']==du.VTYPE_CAP:
            print("rAB{sep}{te}{sep}{v} - z{sep}{te}{sep}{v} <= 0".format(
                sep=du.SEP,
                v=v,
                te=tree_edge
            ),file=out)

def crt35(G,tree_edge,genomes,out):
    for v,data in G.nodes(data=True):
        if data['type']!=du.VTYPE_CAP:
            continue
        if get_genome(G,v)==genomes[0]:
            reps = ["Ab","Aa"]
        else:
            reps = ["Ba","Bb"]
        print("r{a}{sep}{te}{sep}{v} + r{b}{sep}{te}{sep}{v} + y{sep}{te}{sep}{v} >= 1".format(
            a=reps[0],
            b=reps[1],
            te=tree_edge,
            sep=du.SEP,
            v=v
        ),file=out)

def crt36(G,tree_edge,genomes,out):
    for v,data in G.nodes(data=True):
        if data['type']!=du.VTYPE_CAP:
            continue
        if get_genome(G,v)==genomes[0]:
            reps = ["Ab","Aa"]
        else:
            reps = ["Ba","Bb"]
        for r in reps:
            print("y{sep}{te}{sep}{v} + {v} r{r}{sep}{te}{sep}{v} <= {v}".format(
                sep=du.SEP,
                te=tree_edge,
                v=v,
                r=r
            ),file=out)

def crt37(G,tree_edge,genomes,out):
    for x,y,data in G.edges(data=True):
        if data['type'] != du.ETYPE_ADJ:
            continue
        if G.nodes[x]['type'] == du.VTYPE_CAP:
            u,v = y,x
        elif G.nodes[y]['type'] == du.VTYPE_CAP:
            u,v = x,y
        else:
            continue
        if get_genome(G,v) == genomes[0]:
            for r in ["AB","Ab"]:
                print("r{r}{sep}{te}{sep}{v} -  l{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} <= 1".format(
                    r=r,
                    sep=du.SEP,
                    te=tree_edge,
                    v=v,
                    u=u,
                    e=data['id']
                ),file=out)
        else:
            print("rBa{sep}{te}{sep}{v} +  l{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} <= 2".format(
                    sep=du.SEP,
                    te=tree_edge,
                    v=v,
                    u=u,
                    e=data['id']
                ),file=out)

def pcap_reporting(G,tree_edge,genomes,out):
    for c in [crt32,crt33,crt34,crt35,crt36,crt37]:
        c(G,tree_edge,genomes,out)


def csc22(G,tree_edge,out):
    for u,v,data in G.edges(data=True):
        if du.VTYPE_CAP in [G.nodes[x]['type'] for x in [u,v]] or not G.nodes[u].get('cscandidate',False):
            continue
        if data['type']==du.ETYPE_EXTR:
            continue
        print("d{sep}{te}{sep}{v} + d{sep}{te}{sep}{u} + x{sep}{te}{sep}{e} <= 2".format(
            sep=du.SEP,
            te=tree_edge,
            v=v,
            u=u,
            e=data['id']
        ),file=out)
        print("d{sep}{te}{sep}{v} + d{sep}{te}{sep}{u} - x{sep}{te}{sep}{e} >= 0".format(
            sep=du.SEP,
            te=tree_edge,
            v=v,
            u=u,
            e=data['id']
        ),file=out)

def csc23(G,tree_edge,out):
    for u,v,data in G.edges(data=True):
        if data['type']!=du.ETYPE_ID or not G.nodes[u].get('cscandidate',False):
            continue
        print("w{sep}{te}{sep}{u} - w{sep}{te}{sep}{v} = 0".format(
            sep=du.SEP,
            te=tree_edge,
            u=u,
            v=v
        ),file=out)

def csc24(G,tree_edge,out,max_circ_len):
    for u_,v_,data in G.edges(data=True):
        if data['type']!=du.ETYPE_ADJ or du.VTYPE_CAP in [G.nodes[x]['type'] for x in [u_,v_]] or not G.nodes[u_].get('cscandidate',False):
            continue
        for u,v in [(u_,v_),(v_,u_)]:
            print("w{sep}{te}{sep}{u} - w{sep}{te}{sep}{v} + d{sep}{te}{sep}{v} - d{sep}{te}{sep}{u} + {K} x{sep}{te}{sep}{e} - {K} rs{sep}{te}{sep}{u} - {K} rs{sep}{te}{sep}{v} <= {K}".format(
                sep=du.SEP,
                u=u,
                v=v,
                e=data['id'],
                K=max_circ_len,
                te=tree_edge
            ),file=out)

def cs_constraints(G,tree_edge,out,max_circ_len):
    csc22(G,tree_edge,out)
    csc23(G,tree_edge,out)
    csc24(G,tree_edge,out,max_circ_len)



def manual_cs_constraints(circ_singletons,te,out):
    for j, path in enumerate(circ_singletons.values()):
        component_vars = ['x{sep}{te}{sep}{e}'.format(sep=du.SEP,te=te,e=data['id']) for data in path]
        print('{sm} - rms{sep}{te}{sep}{j} <= {totlen}'.format(sm=' + '.join(component_vars), j=j, te=te,sep=du.SEP,
            totlen=len(component_vars)-1), file = out)

def getAllCaps(graphs):
    res = dict((k, set()) for k in set(chain(*graphs.keys())))

    for (child, parent), G in graphs.items():

        for v, vdata in G.nodes(data=True):
            if vdata['type'] == du.VTYPE_CAP:
                res[vdata['id'][0]].add(v)
    return res


# ILP DOMAINS
#

COUNTERS = ['f','n','c','s','pab','pAB','pAb','pAa','pBa','pBb','pABa','pABb']
def domains(graphs, out):
    global_generals = set()
    out.write('bounds\n')
    for te, ((child, parent), G) in enumerate(sorted(graphs.items())):
        for rv in COUNTERS:
            print("0 <= {rv}{sep}{e} <= inf".format(rv=rv,sep=du.SEP,e=te),file=out)
        print("-inf <= q{sep}{e} <= inf".format(sep=du.SEP,e=te))
        for v,data in G.nodes(data=True):
            if data['type']==du.VTYPE_CAP:
                continue
            print("0 <= y{sep}{te}{sep}{v} <= {v}".format(sep=du.SEP,te=te,v=v))
        


def create_sibmap(G,siblings):
    id_to_edges = dict()
    sibmap=dict()
    for u,v,data in G.edges(data=True):
        if data['type']==du.ETYPE_EXTR:
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


            

def warm_start_decomposition(graphs,siblings):
    create_max_weight(graphs)
    for tree_edge, ((child, parent), G) in enumerate(sorted(graphs.items())):
        sibmap=create_sibmap(G,siblings[(child,parent)])
        #cycle_greedy_partial_matching(G,sibmap)
        heuristic_partial_matching(G,sibmap)
        greedy_matching_completion(G,sibmap)


def create_max_weight(graphs):
    root_candidates = dict()
    not_root = set()
    matchings = dict()
    
    for tree_edge, ((child, parent), G) in enumerate(sorted(graphs.items())):
        print(parent,child, file=stderr)
        root_candidates[parent]=child
        print(root_candidates,file=stderr)
        not_root.add(child)
        print(root_candidates,file=stderr)
        matchings[child]=get_max_match(child, G)
    print(root_candidates,file=stderr)
    for x in not_root:
        root_candidates.pop(x,None)
    assert(len(root_candidates)==1)
    parent,child = list(root_candidates.items())[0]
    G=graphs[(child,parent)]
    matchings[parent]=get_max_match(parent,G)
    for tree_edge, ((child, parent), G) in enumerate(sorted(graphs.items())):
        set_based_on_matching(G,matchings[child],child)
        set_based_on_matching(G,matchings[parent],parent)
        #print(G.edges(data=True),file=stderr)
    

def get_anc_map(g):
    return dict(((anc,v) for v,anc in g.nodes(data='anc')))


def set_based_on_matching(G,mtch,genome):
    local = extract_local_graph(genome, G)
    ancmap = get_anc_map(local)
    for u_,v_ in mtch:
        u,v = ancmap[u_],ancmap[v_]
        candidates = [k for k in G[u][v] if G[u][v][k]['type']==du.ETYPE_ADJ]
        assert(len(candidates)==1)
        G[u][v][candidates[0]]['is_set']=True
    for u,v,k,data in local.edges(keys=True,data=True):
        if data['type']!=du.ETYPE_ADJ:
            #skip indel edges
            continue
        if not 'is_set' in G[u][v][k]:
            G[u][v][k]['is_set']=False


def get_max_match(genome, G):
    medges = []
    local = extract_local_graph(genome, G)
    for u,v,tp in list(local.edges(data='type')):
        if tp!=du.ETYPE_ADJ:
            local.remove_edge(u,v)
    telomeres = [u for u,data in G.nodes(data=True) if data['type']==du.VTYPE_CAP]
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
    print("Calculating maximum weight matching",file=stderr)
    mwmtch = list(nx.max_weight_matching(local,weight='mod_weight'))
    print("Maximum weight matching done",file=stderr)
    #print(mwmtch,file=stderr)
    #check matching
    matched = set([v for v,_ in mwmtch]+[v for _,v in mwmtch])
    for v,tp in local.nodes(data='type'):
        if tp == du.VTYPE_CAP:
            continue
        assert(v in matched)
    #mwmtch=greedy_extend_max_match(local,mwmtch)
    #print(mwmtch,file=stderr)
    return [(G.nodes[x]['anc'],G.nodes[y]['anc']) for x,y in mwmtch]


def greedy_extend_max_match(local,mwmtch):
    mwmap = dict()
    for x,y in mwmtch:
        mwmap[x]=y
        mwmap[y]=x
    for v,data in local.nodes(data=True):
        if data['type']==du.VTYPE_CAP:
            #caps don't need to be matched
            continue
        if v in mwmap:
            continue
        #print("Not in Map",file=stderr)
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



def greedy_matching_completion(G,sibmap):
    for v in G.nodes():
        extren = [(v,u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==du.ETYPE_EXTR]
        iden = [(v,u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==du.ETYPE_ID]
        set_extren = [(x,y,z) for x,y,z in extren if G[x][y][z].get('is_set',False)]
        set_iden = [(x,y,z) for x,y,z in iden if G[x][y][z].get('is_set',False)]
        assert(len(set_extren)<=1)
        assert(len(set_iden)<=1)
        assert(len(set_iden)+len(set_extren)<=1)
        if G.nodes[v]['type']==du.VTYPE_CAP:
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
            k__candidates = [k_ for k_ in G[v_][u_] if G[u_][v_][k_]['type']==du.ETYPE_EXTR]
            assert(len(k__candidates)==1)
            k_ = k__candidates.pop()
            corr_at = {u:v,v:u,v_:u_,u_:v_}
            wrongs = set()
            for x,y in corr_at.items():
                wrong = set([(min(w,x),max(w,x),k) for w in G[x] for k in G[x][w] if G[x][w][k]['type']==du.ETYPE_EXTR and w!=y])
                wrongs.update(wrong)
            for x,y,z in wrongs:
                assert(not G[x][y][z].get('is_set',False))
                G[x][y][z]['is_set']=False
            G[u][v][k]['is_set']=True
            G[u_][v_][k_]['is_set']=True



def fill_matching_greedy(G,sibmap):
    for v in G.nodes():
        exten = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==du.ETYPE_EXTR and G[v][u][k].get('is_set',False)]
        assert(len(exten)<=1)
        adjen = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==du.ETYPE_ADJ and G[v][u][k].get('is_set',False)]
        
        iden = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==du.ETYPE_ID and G[v][u][k].get('is_set',False)]
        
        assert(len(adjen)==1)
        assert(len(iden)<=1)
        if G.nodes[v]['type']==du.VTYPE_CAP and len(adjen) == 0:
            #non-active telomere, skip
            continue
        if len(exten) == 0 and len(iden)==0:
            #not matched yet
            if G.nodes[v]['type']==du.VTYPE_CAP:
                continue
            success=False
            for u in G[v]:
                for k in G[v][u]:
                    xt_cands_u = [(x,k) for x in G[u] for k in G[u][x] if G[u][x][k]['type']==du.ETYPE_EXTR and G[u][x][k].get('is_set',False)]
                    if G[v][u][k]['type']==du.ETYPE_EXTR and G[v][u][k].get('is_set',True) and len(xt_cands_u)==0:
                        G[v][u][k]['is_set']=True
                        v_,u_ = sibmap[(v,u)]
                        candidates = [k for k in G[v_][u_]]
                        #extremity edges are always unique
                        assert(len(candidates)==1)
                        G[v_][u_][candidates[0]]['is_set']=True
                        success=True
                        break
            if not success:
                iden = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==du.ETYPE_ID]
                try:
                    assert(len(iden)==1)
                except AssertionError:
                    check_connectedness(G,v)
                    print(v,file=stderr)
                    raise AssertionError("Well ok, at least it's my own problem.")
                #set the indel edge
                u,k=iden[0]
                G[v][u][k]['is_set']=True
        

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
                    if G[v][u][k]['type']==du.ETYPE_EXTR and u not in gnmb:
                        gnmb.add(u)
                        changed=True
        for v in gnmb:
            for u in G[v]:
                for k in G[v][u]:
                    if G[v][u][k]['type']==du.ETYPE_EXTR and u not in gnma:
                        gnma.add(u)
                        changed=True
    assert(len(gnma)<=len(gnmb))
    for x in gnma:
        for y in gnmb:
            xte_candidates = [k for k in G[x][y] if G[x][y][k]['type']==du.ETYPE_EXTR]
            assert(len(xte_candidates)==1)




def cycle_greedy_partial_matching(G,sibmap,max_depth=None):
    trees = dict()
    frontiers = dict()
    depths = dict()
    adjacencies = set(((u,v) for u,v,etype in G.edges(data='type') if etype == du.ETYPE_ADJ and G[u][v].get('is_set',True)))
    for v in G.nodes():
        tree = {du.ETYPE_ADJ : dict(), du.ETYPE_EXTR : dict()}
        frontier = [v]
        trees[v]=tree
        frontiers[v]=frontier
        depths[v]=dict()
    changed = True
    etype = du.ETYPE_ADJ
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
            cand = [k_ for k_ in G[u_][v_] if G[u_][v_][k_]['type']==du.ETYPE_EXTR]
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
                wrong = set([(min(w,x),max(w,x),k) for w in G[x] for k in G[x][w][k] if G[x][w][k]['type']==du.ETYPE_EXTR and w!=y])
                to_remove.update(wrong)
        if len(to_set.intersection(to_remove))>0:
            #failed cycle
            continue
        violated = [(x,y,z) for x,y,z in to_remove if G[x][y][z].get('is_set',False)]
        if len(violated) > 0:
            continue
        for u,v,k in to_set:
            G[u][v][k]['is_set']=True
            




    
def trace_cycle(G,u,v,x,trees,etype):
    xedges = []
    aedges = []
    ae,xe=traceback_single_vertex(G, u, trees,x, etype)
    aedges.extend(ae)
    xedges.extend(xe)
    ae,xe=traceback_single_vertex(G, v, trees,x, invert_etype(etype))
    return aedges,xedges

def traceback_single_vertex(G, u, trees, x, etype):
    xedges=[]
    aedges=[]
    curr = x
    curr_etype = etype
    while curr != u:
        next = trees[u][curr_etype][curr]
        candidates = [(curr,next,k) for k in G[curr][next] if G[curr][next][k]==curr_etype]
        assert(len(candidates)==1)
        if curr_etype==du.ETYPE_ADJ:
            aedges.extend(candidates)
        else:
            xedges.extend(candidates)
    return aedges,xedges


def heuristic_partial_matching(G,sibmap,max_iter=None):
    #TODO: Implement max_iter
    trees = dict()
    frontiers = dict()
    #TODO: Filter out useless adjacencies, i.e. those that don't have extremity-edges
    adjacencies = set(((u,v) for u,v,etype in G.edges(data='type') if etype == du.ETYPE_ADJ and G[u][v].get('is_set',True)))
    for v in G.nodes():
        tree = {du.ETYPE_ADJ : dict(), du.ETYPE_EXTR : dict()}
        frontier = [v]
        trees[v]=tree
        frontiers[v]=frontier
    #bfs for each node starting with extremity edge until end of a path is reached or connection to neighboring adjacency found
    changed = True
    etype = du.ETYPE_ADJ
    
    while changed:
        etype = invert_etype(etype)
        rm_adj = set()
        changed = False
        print('-----',file=stderr)
        for v in G.nodes():
            extren = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==du.ETYPE_EXTR and G[v][u][k].get('is_set',False)]
            print(v,extren,file=stderr)
            assert(len(extren)<=1)
        for u,v in adjacencies:
            bad=False
            bfs = bfs_step(G, trees, frontiers, etype, u)
            changed = changed or bfs
            bfs = bfs_step(G, trees, frontiers, etype, v)
            changed = changed or bfs
            if len(trees[u][du.ETYPE_EXTR])==0 or len(trees[v][du.ETYPE_EXTR]) == 0:
                rm_adj.add((u,v))
                bad=True
            else:
                #TODO: comment out
                for x in (u,v):
                    print(x,trees[x],file=stderr)
                    assert(x in trees[x][du.ETYPE_EXTR].values())
            #if bf != changed:
            #    print(frontiers[u],file=stderr)
            #    print(trees[v][invert_etype(etype)],file=stderr)
            #    print("-------------",file=stderr)
            #check whether anything special happened at the frontiers, i.e. telomere, indel, cycle
            if not bad:
                for x in frontiers[u]:
                    if x in trees[v][invert_etype(etype)]:
                        #cycle detected
                        for x in (u,v):
                            assert(x in trees[x][du.ETYPE_EXTR].values())
                        x,aedges=try_fix_cycle(G,sibmap,u,v,trees,frontiers,x,etype)
                        if x:
                            #remove the adjacency edges now dealt with
                            rm_adj.update(aedges)
                            break
                        for x in (u,v):
                            assert(x in trees[x][du.ETYPE_EXTR].values() or len(trees[x][du.ETYPE_EXTR])==0)
        adjacencies.difference_update(rm_adj)
    for v in G.nodes():
            extren = [(u,k) for u in G[v] for k in G[v][u] if G[v][u][k]['type']==du.ETYPE_EXTR and G[v][u][k].get('is_set',False)]
            print(v,extren,file=stderr)
            assert(len(extren)<=1)

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
        if curr_etype==du.ETYPE_EXTR:
            candidates = [k for k in G[curr][x_] if G[curr][x_][k]['type']==du.ETYPE_EXTR]
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
        print("Tf cycle loop B {} curr {}".format(count,curr),file=stderr)
        print(v,file=stderr)
        print(trees[v],file=stderr)
        try:
            x_=trees[v][curr_etype][curr]
        except KeyError:
            success=False
            break
        if curr_etype==du.ETYPE_EXTR:
            candidates = [k for k in G[curr][x_] if G[curr][x_][k]['type']==du.ETYPE_EXTR]
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
    

def try_fix_edges(G,edges,sibmap):
    to_set = set()
    to_remove = set()
    all_edges=set()
    for x,y in edges:
        to_set.add(tuple(sorted((x,y))))
        x_,y_ = sibmap[(x,y)]
        to_set.add(tuple(sorted((x_,y_))))
        #Remove other edges
        rm = [tuple(sorted((w,z))) for w in [x,y,x_,y_] for z in G[w] for k in G[w][z] if G[w][z][k]['type']==du.ETYPE_EXTR and z not in [x,y,x_,y_]]
        rm_sib = [tuple(sorted(sibmap[e])) for e in rm]
        to_remove.update(rm)
        to_remove.update(rm_sib)
        all_edges.update([tuple(sorted((w,z))) for w in [x,y,x_,y_] for z in G[w] for k in G[w][z] if G[w][z][k]['type']==du.ETYPE_EXTR])
    try:
        assert(len(to_remove.union(to_set).symmetric_difference(all_edges))==0)
    except AssertionError:
        print("To remove: ",to_remove,file=stderr)
        print("To set: ",to_set,file=stderr)
        print("All: ",all_edges,file=stderr)
        raise AssertionError("Did not consider all edges for all vertices :/")
    if len(to_set.intersection(to_remove))>0:
        #bad, conflicting edges :/
        return False
    for x,y in to_remove:
        candidates = [k for k in G[x][y] if G[x][y][k]['type']==du.ETYPE_EXTR]
        assert(len(candidates)==1)
        if  G[x][y][candidates[0]].get('is_set',False):
            return False
    for x,y in to_set:
        candidates = [k for k in G[x][y] if G[x][y][k]['type']==du.ETYPE_EXTR]
        assert(len(candidates)==1)
        G[x][y][candidates[0]]['is_set']=True
    for x,y in to_remove:
        candidates = [k for k in G[x][y] if G[x][y][k]['type']==du.ETYPE_EXTR]
        assert(len(candidates)==1)
        G[x][y][candidates[0]]['is_set']=False
        #print("Removed" , x,y,"because of",to_set,file=stderr)
    

def remove_from_tree(trees,u,x,x_,etype):
    curr = x
    curr_etype=etype
    while curr != x_:
        curr=trees[u][curr_etype].pop(curr)
        curr_etype=invert_etype(curr_etype)


def cleanup_tree(G,trees,frontiers,u):
    return
    print("Change no",u,"cleanup",trees[u],file=stderr)
    #convert to a parent -> children tree
    pc = {du.ETYPE_ADJ:dict(),du.ETYPE_EXTR:dict()}
    if len(trees[u][du.ETYPE_EXTR])==0:
        #nothing to be done, reset adjacencies and frontiers
        trees[u][du.ETYPE_ADJ]=dict()
        frontiers[u]=[]
        return
    for etype in [du.ETYPE_ADJ,du.ETYPE_EXTR]:
        for child,parent in trees[u][etype].items():
            if not parent in pc[etype]:
                pc[etype][parent]=[]
            pc[etype][parent].append(child)
    print('ä'*20,file=stderr)
    print(u,trees[u],file=stderr)
    print(u,pc,file=stderr)
    print(u,trees[u][du.ETYPE_EXTR].values(),file=stderr)
    print('ä'*20,file=stderr)
    assert(u in trees[u][du.ETYPE_EXTR].values())
    assert(u in pc[du.ETYPE_EXTR])
    stack = [(u,du.ETYPE_EXTR)]
    while len(stack)>0:
        x,etype = stack.pop()
        for child in pc[etype].get(x,[]):
            candidates = [k for k in G[x][child] if G[x][child][k]['type']==etype and G[x][child][k].get('is_set',True)]
            if len(candidates) == 0:
                #dead edge, remove
                remove_dead_branch(trees,frontiers,u,pc,child,invert_etype(etype))
            else:
                stack.append((child,invert_etype(etype)))
    print("Change no ",u,"cleanup done",trees[u],file=stderr)



def remove_dead_branch(trees,frontiers,u,pc,start,petype):
    dead = {du.ETYPE_ADJ : [], du.ETYPE_EXTR : []}
    stack = [(start,petype)]
    print('*'*20,file=stderr)
    print("Before",file=stderr)
    print(start,petype,file=stderr)
    print(trees[u],file=stderr)
    print(pc,file=stderr)
    print('*'*20,file=stderr)
    while len(stack)>0:
        x,etype = stack.pop()
        #print("Stacklen", len(stack),file=stderr)
        dead[invert_etype(etype)].append(x)
        for child in pc[etype].get(x,[]):
            stack.append((child,invert_etype(etype)))
    for etype, deads in dead.items():
        frontiers_new = set(frontiers[u]).difference(set(deads))
        frontiers[u]=frontiers_new
        for d in deads:
            print(u,etype,d,file=stderr)
            trees[u][etype].pop(d)
    print('*'*20,file=stderr)
    print("After",file=stderr)
    print(start,petype,file=stderr)
    print(trees[u],file=stderr)
    print(pc,file=stderr)
    print('*'*20,file=stderr)




def invert_etype(etype):
    return du.ETYPE_EXTR if etype==du.ETYPE_ADJ else du.ETYPE_ADJ
            

def bfs_step(G, trees, frontiers, etype, u,forbidden=None):
    if forbidden is None:
        forbidden=[u]
    newfrontieru = []
    changed = False
    #print(etype,file=stderr)
    for x in frontiers[u]:
        for x_ in G[x]:
            for k in G[x][x_]:
                if G[x][x_][k]['type']==etype and G[x][x_][k].get('is_set',True) and x_ not in trees[u][du.ETYPE_ADJ] and x_ not in trees[u][du.ETYPE_EXTR] and x_ not in forbidden:
                    print("Found new edge! {} -> {}".format(x,x_),file=stderr)
                    print("Change",u,trees[u],file=stderr)
                    newfrontieru.append(x_)
                    trees[u][etype][x_]=x
                    changed=True
    oldfrontiers = frontiers[u]
    frontiers[u]=newfrontieru
    if changed:
        assert(u in trees[u][du.ETYPE_EXTR].values())
    return changed


        
#COUNTERS = ['f','n','c','s','pab','pAB','pAb','pAa','pBa','pBb','pABa','pABb']


def sol_from_decomposition(graphs,alpha,out):
    #TODO: Implement
    var_map = dict()
    wsum = 0.0
    fsum = 0
    for tree_edge, ((child, parent), G) in enumerate(sorted(graphs.items())):
        genomes = [child,parent]
        local = G.copy()
        print('\n'.join([str(d) for d in local.nodes(data=True)]),file=stderr)
        print("---------****___",file=stderr)
        print('\n'.join([str(d) for d in local.edges(data=True)]),file=stderr)
        counts = {'2n':0,'w':0.0,'c':0,'ab':0,'AB':0,'Aa':0,'Ab':0,'Ba':0,'Bb':0}
        for u,v,k,d in G.edges(keys=True,data=True):
            xvar = "x{sep}{te}{sep}{e}".format(sep=du.SEP,te=tree_edge,e=d['id'])
            if not d.get('is_set',False) or d['type']==du.ETYPE_ID:
                local.remove_edge(u,v,key=k)
            var_map[xvar]= 1 if d.get('is_set',False) else 0
            if d['type']==du.ETYPE_ADJ:
                ancvar = 'a{sep}{e}'.format(sep=du.SEP,e=G[u][v][k]['anc'])
                if not ancvar in var_map:
                    var_map[ancvar] =  1 if d.get('is_set',False) else 0
                else:
                    assert(var_map[ancvar] ==  (1 if d.get('is_set',False) else 0))
                if d.get('is_set',False):
                    counts['w']+=d['weight']
            if d['type']==du.ETYPE_EXTR and d.get('is_set',False): 
                counts['2n']+=1
            assert(local.has_edge(u,v,key=k)==(var_map[xvar]==1) or d['type']==du.ETYPE_ID)
        print('*'*20,file=stderr)
        print('\n'.join([str(d) for d in local.edges(data=True)]),file=stderr)
        print('*'*20,file=stderr)
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
                print(comp)
                v=list(comp)[0]
                print(v,G.nodes[v]['type'])
                assert(G.nodes[path_ends[0]]['type']==du.VTYPE_CAP)
                continue
            assert(len(path_ends)==2 or len(path_ends)==0)
            is_positive=True
            print("path ends",tree_edge,path_ends,file=stderr)
            for v in path_ends:
                if G.nodes[v]['type']!=du.VTYPE_CAP:
                    is_positive=False
            if is_positive:
                for v in comp:
                    var_map["y{sep}{te}{sep}{v}".format(v=v,sep=du.SEP,te=tree_edge)]=min_id
                    var_map["z{sep}{te}{sep}{v}".format(v=v,sep=du.SEP,te=tree_edge)]= 0 if v!=min_id else 1
            else:
                for v in comp:
                    var_map["y{sep}{te}{sep}{v}".format(v=v,sep=du.SEP,te=tree_edge)]=0
                    var_map["z{sep}{te}{sep}{v}".format(v=v,sep=du.SEP,te=tree_edge)]=0
            if len(path_ends)==0:
                #cycle
                counts['c']+=1
                for v in comp:
                    var_map["l{sep}{te}{sep}{v}".format(v=v,sep=du.SEP,te=tree_edge)]= 0
                    if get_genome(G,v) == genomes[0]:
                        var_map["rc{sep}{te}{sep}{v}".format(sep=du.SEP,te=tree_edge,v=v)] = 0 if v!=min_id else 1
                        var_map["rab{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge)]=0
            else:
                mn = min(path_ends)
                mx = max(path_ends)
                if get_genome(G,mn) == get_genome(G,mx):
                #same end, so can just pick this as the genome
                    for v in comp:
                        var_map["l{sep}{te}{sep}{v}".format(v=v,sep=du.SEP,te=tree_edge)]= 0 if get_genome(G,mx) == genomes[0] else 1
                else:
                    mxg = 0 if get_genome(G,mx)==genomes[0] else 1
                    mng = 0 if get_genome(G,mn)==genomes[0] else 1
                    for v in comp:    
                        var_map["l{sep}{te}{sep}{v}".format(v=v,sep=du.SEP,te=tree_edge)]= mxg if v!=mn else mng
                print(comp,file=stderr)
                for v in comp:
                    if v in [mx,mn] or get_genome(G,v)!=genomes[0]:
                        continue
                    var_map["rab{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge)]=0
                #set reporting variables at ends
                #TODO: What about even paths?
                if G.nodes[mn]['type']!=du.VTYPE_CAP:
                    #must be type ab, because caps generally have lower ids
                    assert(G.nodes[mx]['type']!=du.VTYPE_CAP)
                    assert(get_genome(G,mn)==genomes[0] or get_genome(G,mn)==get_genome(G,mx))
                    print(path_ends,'potentially pab')
                    if get_genome(G,mn)==genomes[0]:
                        if get_genome(G,mn)!=get_genome(G,mx):
                            print("mn",mn,"mx",mx,file=stderr)
                            var_map["rab{sep}{te}{sep}{v}".format(sep=du.SEP,v=mn,te=tree_edge)]=1
                            counts['ab']+=1
                        else:
                            var_map["rab{sep}{te}{sep}{v}".format(sep=du.SEP,v=mn,te=tree_edge)]=0
                    #do not uncomment this var_map["rab{sep}{te}{sep}{v}".format(sep=du.SEP,v=mx,te=tree_edge)]=0
                else:
                    assert(min_id==mn)
                    if G.nodes[mx]['type']!=du.VTYPE_CAP:
                        #only one var to set
                        if get_genome(G,mx)==genomes[0]:
                            var_map["rab{sep}{te}{sep}{v}".format(sep=du.SEP,v=mx,te=tree_edge)]=0
                    else:
                        if get_genome(G,mx)==genomes[0]:
                            reps = ["Ab","Aa","AB"]
                        else:
                            reps = ["Ba","Bb"]
                        for r in reps:
                            var_map["r{r}{sep}{te}{sep}{v}".format(r=r,sep=du.SEP,v=mx,te=tree_edge)]=0
                    if get_genome(G,mn)==genomes[0]:
                        reps = ["Ab","Aa","AB"]
                    else:
                        reps = ["Ba","Bb"]
                    for r in reps:
                        var_map["r{r}{sep}{te}{sep}{v}".format(r=r,sep=du.SEP,v=mn,te=tree_edge)]=1 if fits_ptype(G,genomes,mx,r[1]) else 0
                        counts[r]+=var_map["r{r}{sep}{te}{sep}{v}".format(r=r,sep=du.SEP,v=mn,te=tree_edge)]
            if len(path_ends)>0:
                mnreportv = ["r{tp}{sep}{te}{sep}{v}".format(sep=du.SEP,te=tree_edge,v=mn,tp=tp) for tp in ['AB','Aa','Ab','Ba','Bb','ab']]
                mnreportvl = [var_map.get(rv,None) for rv in mnreportv]
                mnl = var_map["l{sep}{te}{sep}{v}".format(sep=du.SEP,te=tree_edge,v=mn)]
                mxl = var_map["l{sep}{te}{sep}{v}".format(sep=du.SEP,te=tree_edge,v=mx)]
                print("path ends reported",tree_edge,path_ends,mnreportvl,file=stderr)
                assert(1 in mnreportvl or (mxl == mnl and (G.nodes[mn]['type'],get_genome(G,mn)) == (G.nodes[mx]['type'],get_genome(G,mx))))
        n = int(ceil(counts['2n']/2))
        var_map['n{sep}{te}'.format(sep=du.SEP,te=tree_edge)]=n
        var_map['c{sep}{te}'.format(sep=du.SEP,te=tree_edge)]=counts['c']
        var_map['w{sep}{te}'.format(sep=du.SEP,te=tree_edge)]=counts['w']
        wsum+=counts['w']
        for r in ['ab','AB','Aa','Ab','Ba','Bb']:
            var_map["p{r}{sep}{te}".format(r=r,sep=du.SEP,te=tree_edge)]=counts[r]
        pABa = max(counts['Aa'],counts['Ba'])
        pABb = max(counts['Ab'],counts['Bb'])
        var_map['pABa{sep}{te}'.format(sep=du.SEP,te=tree_edge)]=pABa
        var_map['pABb{sep}{te}'.format(sep=du.SEP,te=tree_edge)]=pABb
        q = int(ceil((counts['ab']+pABa+pABb - counts['AB'])/2))
        var_map['q{sep}{te}'.format(sep=du.SEP,te=tree_edge)]=q
        #TODO: Circular singletons!
        f = n - counts['c'] + q
        var_map['f{sep}{te}'.format(sep=du.SEP,te=tree_edge)]=f
        fsum+=f
    obj = alpha*fsum + (alpha-1)*wsum
    print("# Objective value = {}".format(obj),file=out)
    for vr,vl in var_map.items():
        print(vr,vl,file=out)



def fits_ptype(G,genomes,mx,r):
    if r.isupper() != (G.nodes[mx]['type']==du.VTYPE_CAP):
        return False
    if (r.upper()=='A') != (get_genome(G,mx) == genomes[0]):
        return False
    return True
#
# ILP VARIABLES
#

def variables(graphs,circ_sings, out):

    #
    # integer variables
    #
    global_generals = set()
    out.write('generals\n')
    for te, ((child, parent), G) in enumerate(sorted(graphs.items())):
        LOG.info(('writing general variables for relational diagram of {} ' + \
                'and {}').format(child, parent))
        for rv in COUNTERS:
            print("{rv}{sep}{e}".format(rv=rv,sep=du.SEP,e=te),file=out)
        print("q{sep}{e}".format(sep=du.SEP,e=te))
        for v,data in G.nodes(data=True):
            print("y{sep}{te}{sep}{v}".format(sep=du.SEP,te=te,v=v),file=out)
            if data['type']==du.VTYPE_CAP or not G.nodes[v].get('cscandidate',False):
                continue
            if (child, parent) not in circ_sings:
                print("w{sep}{te}{sep}{v}".format(sep=du.SEP,te=te,v=v),file=out)


    for gg in global_generals:
        print(gg,file=out)
    #
    # binary variables
    #
    out.write('binaries\n')

    global_binaries = set()
    for tree_edge, ((child, parent), G) in enumerate(sorted(graphs.items())):

        LOG.info(('writing binary variables for relational diagram of {} ' + \
                'and {}').format(child, parent))
        genomes = [child,parent]
        if (child, parent) in circ_sings:
            for j in range(len(circ_singletons[(child, parent)])):
                print('rms{sep}{te}{sep}{j}'.format(sep=du.SEP,te=tree_edge,j=j))
        for v,data in G.nodes(data=True):
            print("z{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
            print("l{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
            global_binaries.add("g{sep}{v}".format(sep=du.SEP,v=data['anc']))
            if data['type']!=du.VTYPE_CAP:
                if get_genome(G,v) == genomes[0]:
                    print("rab{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
                    print("rc{sep}{te}{sep}{v}".format(sep=du.SEP,te=tree_edge,v=v),file=out)
                if (child, parent) not in circ_sings:
                    if G.nodes[v].get('cscandidate',False):
                        print("d{sep}{te}{sep}{v}".format(sep=du.SEP,te=tree_edge,v=v),file=out)
                        print("rs{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
            else:
                if get_genome(G,v) == genomes[0]:
                    print("rAB{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
                    print("rAa{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
                    print("rAb{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
                else:
                    print("rBa{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
                    print("rBb{sep}{te}{sep}{v}".format(sep=du.SEP,v=v,te=tree_edge),file=out)
        for u,v,data in G.edges(data=True):
            print("x{sep}{te}{sep}{e}".format(sep=du.SEP,te=tree_edge,e=data['id']),file=out)

    print('\n'.join(global_binaries), file = out)
    print('\n')


def add_weight_tels(G,add_telweight):
    for u,v,k,data in G.edges(keys=True,data=True):
        if data['type']==du.ETYPE_ADJ and du.VTYPE_CAP in [G.nodes[x]['type'] for x in [u,v]]:
            G[u][v][k]['weight']=G[u][v][k]['weight']+add_telweight


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
            G.add_nodes_from(reduce(lambda x, y: x + y, (((g, du.EXTR_HEAD), (g,
                du.EXTR_TAIL)) for g in genes)))
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

#        genes_edg = [((g, du.EXTR_HEAD), (g, du.EXTR_TAIL)) for g in genes]
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
        LOG.info('identified %s candidate telomeres in genome %s' %(
            len(telomeres), species))
    return res

def cmd_arguments(parser):
    parser.add_argument('tree', type=open,
            help='phylogenetic tree as parent-child relation table')
    parser.add_argument('candidateAdjacencies', type=open,
            help='candidate adjacencies of the genomes in the phylogeny')
    parser.add_argument('-t', '--no_telomeres', action='store_true',
            help='don\'t add any additional telomeres. Overrides -at.')
    parser.add_argument('-l', '--all_telomeres', action='store_true',
            help='Add additional telomeres to all nodes.')
    parser.add_argument('-m', '--output_id_mapping', type=FileType('w'),
            help='writs a table with ID-to-gene extremity mapping')
    parser.add_argument('-a', '--alpha', default = 0.5, type=float,
            help='linear weighting factor for adjacency weights vs DCJ ' + \
                    'indel distance (alpha = 1 => maximize only DCJ indel ' + \
                    'distance)')
    parser.add_argument('-b', '--beta', type=float,
            help='Backwards compatible beta parameter from v1. Telomeres will be re-weighted and the ILP scaled to be equivalent to v1.')
    parser.add_argument('--def-telomere-weight', '-w', default = 0, type=float,
            help='Default weight for added telomeres. Has no effect if -t is used. For most use cases this should be <= 0.')
    parser.add_argument('--set-circ-sing-handling',choices=['adaptive','enumerate','barbershop'],default='adaptive')
    parser.add_argument('-s', '--separator', default = du.DEFAULT_GENE_FAM_SEP, \
            help='Separator of in gene names to split <family ID> and ' +
                    '<uniquifying identifier> in adjacencies file')
    parser.add_argument('-ws','--warm-start-sol')
    parser.add_argument('-plb','--pairwise-lower-bnds')

def main(args):
    # load & process input data
    LOG.info('loading species tree from {}'.format(args.tree.name))
    speciesTree = du.parseTree(args.tree)

    LOG.info(('loading candidate adjacencies from {}, using "{}" to separate' + \
            ' gene family from uniquifying identifier').format(
        args.candidateAdjacencies.name, args.separator))
    candidateAdjacencies = du.parseAdjacencies(args.candidateAdjacencies,
                                               sep=args.separator)
    leaves = set([x for x, v in du.getLeaves(speciesTree).items() if v])
    LOG.info("Leaves of phylogeny: {}".format(leaves))
    # add telomeres
    telomeres = identifyCandidateTelomeres(candidateAdjacencies,args.def_telomere_weight,leaves, dont_add=args.no_telomeres,addToAll=args.all_telomeres)

    # construct adjacency graphs
    genes = candidateAdjacencies['genes']
    adjacencies = candidateAdjacencies['adjacencies']
    weights = candidateAdjacencies['weights']

    global_ext2id = du.IdManager(0,is_cap=lambda x: False)
    LOG.info(('constructing relational diagrams for all {} branches of ' + \
            'the tree').format(len(speciesTree)))
    relationalDiagrams = du.constructRelationalDiagrams(speciesTree,
            adjacencies, telomeres, weights, genes, global_ext2id,
            sep=args.separator)

    graphs = relationalDiagrams['graphs']

    siblings = relationalDiagrams['siblings']

    for G in graphs.values():
        du.checkGraph(G,cf=True,checkForAllTels=args.all_telomeres and not args.no_telomeres)
        #if 7 in G:
        #    LOG.info(G[7])

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
            LOG.info(f'identified {len(circ_singletons[ident])} circular singleton candidates')
        elif max_tolerable > 0:
            try:
                circ_singletons[ident] = du.identifyCircularSingletonCandidates(G,max_number=max_tolerable)
                LOG.info(f'identified {len(circ_singletons[ident])} circular singleton candidates')
            except du.TooManyCSException:
                pass


    #caps = getAllCaps(graphs)
    # construct & output ILP
    out = stdout
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
        LOG.info('Re-weighting telomeres: s={s},add_tel_weight={adt},recalculated_alpha={recalph}'.format(s=scale_factor,adt=our_beta_add_telweight,recalph=our_alpha))
        for ident,G in graphs.items():
            add_weight_tels(G,our_beta_add_telweight)
    lbounds = {}
    if args.pairwise_lower_bnds:
        LOG.info("Reading Lower bound matrix")
        lbounds = du.parse_lower_bound_file(args.pairwise_lower_bnds)


    LOG.info('writing objective over all graphs')
    objective(graphs, our_alpha, out)

    LOG.info('writing constraints...')
    constraints(graphs, siblings,circ_singletons, out,lower_bound_mat=lbounds)

    LOG.info('writing domains...')
    domains(graphs, out)

    LOG.info('writing variables...')
    variables(graphs,circ_singletons,out)

    if args.output_id_mapping:
        LOG.info('writing ID-to-gene extremity mapping to {}'.format(
            args.output_id_mapping.name))
        idMap = global_ext2id.getMap()
        out_table = list()
        for k, v in idMap.items():
            out_table.append((str(v), k[0], k[1][0], k[1][1]))
        out_table.sort(key = lambda x: int(x[0]))
        print('\n'.join(map(lambda x: '\t'.join(x), out_table)),
                file=args.output_id_mapping)

    LOG.info('DONE')
    out.write('end\n')
    if args.warm_start_sol:
        LOG.info('run warm start')
        warm_start_decomposition(graphs,siblings)
        with open(args.warm_start_sol,'w') as f:
            sol_from_decomposition(graphs,our_alpha,f)
        LOG.info('DONE')
    

