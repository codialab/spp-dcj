# import from built-in packages
from argparse import FileType
from sys import stdout
from itertools import chain
from functools import reduce
import logging

# import from third-party packages

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

def constraints(graphs,families,fam_bounds, siblings,circ_singletons, out,lower_bound_mat={}):
    out.write('subject to\n')
    for i, ((child, parent), G) in enumerate(sorted(graphs.items())):

        LOG.info(('writing constraints for relational diagram of {} and ' + \
                '{}').format(child, parent))
        genomes = [child,parent]
        global_constraints(G,families,fam_bounds,out)
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


def c01(G,families,fam_bounds,out):
    id_to_v = dict()
    ht_pairs = []
    for v,data in G.nodes(data=True):
        if data['type']==du.VTYPE_CAP:
            continue
        vid = data['id']
        ovid = du.complement_id(vid)
        if ovid in id_to_v:
            ht_pairs.append((v,id_to_v.pop(ovid)))
        else:
            id_to_v[vid]=v
    assert(len(id_to_v)==0)
    for u,v in ht_pairs:
        f=du.nodeGetFamily(G,u)
        assert(f==du.nodeGetFamily(G,v))
        gnm = get_genome(G,v)
        if len(families[gnm][f])==fam_bounds[gnm][f][0]:
            #lower bound is the same as the family size, everything has to be set
            print("g{sep}{v} = 1".format(sep=du.SEP,v=G.nodes[v]['anc']),file=out)
            print("g{sep}{u} = 1".format(sep=du.SEP,u=G.nodes[u]['anc']),file=out)
        else:
            print("g{sep}{v} -  g{sep}{u} = 0".format(sep=du.SEP,v=G.nodes[v]['anc'],u=G.nodes[u]['anc']),file=out)


def c02_new(G,fam_bounds,out):
    family_heads = dict()
    for v, data in G.nodes(data=True):
        if data['type']==du.VTYPE_CAP:
            continue
        vid = data['id']
        if vid[1][1] == du.EXTR_HEAD:
            gnm = du.get_genome(G,v)
            if not gnm in family_heads:
                family_heads[gnm] = dict()
            fam = du.nodeGetFamily(G,v)
            if not fam in family_heads[gnm]:
                family_heads[gnm][fam]=[]
            family_heads[gnm][fam].append(data['anc'])
    for g, fhds in family_heads.items():
        for f,hds in fhds.items():
            sm = " + ".join(["g{sep}{v}".format(sep=du.SEP,v=v) for v in hds])
            print("{sm} <= {hb}".format(sm=sm, hb=fam_bounds[g][f][1]),file=out)
            print("{sm} >= {lb}".format(sm=sm, lb=fam_bounds[g][f][0]),file=out)
    
    #print("g{sep}{v} = 1".format(sep=du.SEP,v=data['anc']),file=out)

#Old,removed
#def c02(G,out):
#    for u,v in get_gene_extremities(G):
#        print("g{sep}{u} - g{sep}{v} = 0".format(u=G.nodes[u]['anc'],v=G.nodes[v]['anc'],sep=du.SEP),file=out)

def c02(G,out):
    for v,data in G.nodes(data=True):
        aes = [G[v][u][i]['anc'] for u in G[v] for i in G[v][u] if G[v][u][i]['type']==du.ETYPE_ADJ]
        sm = ' + '.join(['a{sep}{e}'.format(sep=du.SEP,e=e) for e in aes])
        print("{sm} - g{sep}{v} = 0".format(sm=sm,sep=du.SEP,v=G.nodes[v]['anc']),file=out)


def global_constraints(G,families,fam_bounds,out):
    c01(G,families,fam_bounds,out)
    c02(G,out)

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
    c23_new(G,genomes,tree_edge,out)

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


def check_b_necessary(G,u):
    extr_cands = [(x,k) for x in G[u] for k in G[u][x] if G[u][x][k]['type']==du.ETYPE_EXTR]
    if len(extr_cands)==0:
            #other genome does not have this gene, mm is automatic
            return False
    x,_ = extr_cands.pop()
    ide_cands = [(y,k) for y in G[x] for k in G[x][y] if G[x][y][k]['type']==du.ETYPE_ID]
    if len(ide_cands)==0:
            #other genome cannot delete, mm is automatic
            return False
    return True

def c23_new(G,genomes,tree_edge,out):
    for u,v,data in G.edges(data=True):
        if data['type']!=du.ETYPE_ID:
            continue
        f = du.nodeGetFamily(G,u)
        #check if we even have to set this
        if not check_b_necessary(G,u):
            continue
        if du.get_genome(G,u)==genomes[0]:
            print("x{sep}{te}{sep}{e} - b{sep}{te}{sep}{f} <= 0".format(sep=du.SEP,te=tree_edge,e=data['id'],f=f),file=out)
        else:
            print("x{sep}{te}{sep}{e} + b{sep}{te}{sep}{f} <= 1".format(sep=du.SEP,te=tree_edge,e=data['id'],f=f),file=out)

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
        bmvars = set()
        for u,v,data in G.edges(data=True):
            print("x{sep}{te}{sep}{e}".format(sep=du.SEP,te=tree_edge,e=data['id']),file=out)
            if data['type']==du.ETYPE_ID:
                if check_b_necessary(G,u):
                    f=du.nodeGetFamily(G,u)
                    bmvars.add("b{sep}{te}{sep}{f}".format(sep=du.SEP,te=tree_edge,f=f))
        for b in bmvars:
            print(b,file=out)

    print('\n'.join(global_binaries), file = out)
    print('\n')


def cmd_arguments(parser):
    parser.add_argument('tree', type=open,
            help='phylogenetic tree as child->parent relation table (Important: children must be in first column.)')
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
    parser.add_argument('--family-bounds','-fmb',type=open)

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

    global_ext2id = du.IdManager(0,is_cap=lambda x: False)
    LOG.info(('constructing relational diagrams for all {} branches of ' + \
            'the tree').format(len(speciesTree)))
    relationalDiagrams = du.constructRelationalDiagrams(speciesTree,
            adjacencies, telomeres, weights, genes, global_ext2id,fam_bounds=fam_bounds,
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
    constraints(graphs, families,fam_bounds,siblings,circ_singletons, out,lower_bound_mat=lbounds)

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
        for (child,parent) in relationalDiagrams['idmanagers']:
            loc_idmap = relationalDiagrams['idmanagers'][(child,parent)].getMap()
            for k,v in loc_idmap.items():
                print('\t'.join([child,parent,str(v),k[0],k[1][0],k[1][1]]),file=args.output_id_mapping)
    LOG.info('DONE')
    out.write('end\n')
    if args.warm_start_sol:
        LOG.info('run warm start')
        du.warm_start_decomposition(graphs,fam_bounds,families,siblings)
        with open(args.warm_start_sol,'w') as f:
            du.sol_from_decomposition(graphs,circ_singletons,our_alpha,f)
        LOG.info('DONE')
    

