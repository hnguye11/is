from __future__ import division
import networkx as nx
from copy import copy
import random
from random import random as rnd, randint as rndi
from math import log, exp, sqrt
import matplotlib.pyplot as plt
import numpy as np
import time

from util import *


def PRINT_SIM_INFO_DETAILED(ns, i, mu_h, pp, x):
    if SIM_INFO != "": print SIM_INFO
    print "ns, i  = %d, %d"%(ns, i)
    print "mu_h   = %e, %e"%(mu_h[0], mu_h[1])
    print "pi, pp = %e, %e"%(p[i], pp)
    print "x[%02d]  = %d"%(i, x[i])
    print "cache  = %d (hit:%d)"%(len(CSET) - 2, CSET["hit"])
    print "x      = %s"%(TO_STRING(x[:i+1]) + "_"*(M-i-1))
    print
    
    
def PRINT_SIM_INFO(ns, x, w):
    if SIM_INFO != "": print SIM_INFO
    print "ns     =", ns
    print "L(x)   =", CALC_LOSS(x)
    print "w(x)   = %e"%w
    print "cache  = %d (hit:%d)"%(len(CSET) - 2, CSET["hit"])
    print "x      =", TO_STRING(x)
    print
    

def CALC_LOSS(x):
    G = nx.DiGraph()
    G.add_nodes_from(V)
    for i,Ei in enumerate(E):
        if x[i] == 1: G.add_edge(Ei[0], Ei[1])

    Vx = [V[0]] + list(nx.descendants(G, V[0]))
    return sum(L[n] for n in range(N) if V[n] in Vx)


def VERTEX_SEARCH(x, M0):    
    Gx = nx.DiGraph()
    
    for i in range(M0):
        if x[i] == 1: Gx.add_edge(E[i][0], E[i][1], weight=0.0)
        
    # Vx1 contains all vertices that are guaranteed reachable from s
    Vx1 = [] if not V[0] in Gx.nodes() else list(nx.descendants(Gx, V[0]))
    
    for i in range(M0,M):
        Gx.add_edge(E[i][0], E[i][1], weight=logp[i])
    assert(V[0] in Gx.nodes())

    # Vx2 contains all vertices that could be reachable from s
    Vx2 = list(nx.descendants(Gx, V[0]))
    assert(all(Vi in Vx2 for Vi in Vx1))
    
    # removing nodes in Gx that are not in Vx2
    for Vi in Gx.nodes():
        if Vi != V[0] and not Vi in Vx2: Gx.remove_node(Vi)
    
    # initial setup for vertex search
    Vx = [Vi for Vi in Vx2 if not Vi in Vx1]
    Nx = len(Vx)
    SearchTree = [([1] * Nx, 0)]
    px = 0.0
    t = 0
    
    while t < len(SearchTree):
        y, j = SearchTree[t]
        
        if j < Nx:
            SearchTree.append((y, j+1))
            z = y[:j] + [0] + [1] * (Nx-j-1)
            Vz = [V[0]] + Vx1 + [Vx[k] for k in range(Nx) if z[k]==1]
            
            if sum(L[n] for n in range(N) if V[n] in Vz) > L_THRES:
                SearchTree.append((z, j+1))

        elif j == Nx:
            Vy = [V[0]] + Vx1 + [Vx[k] for k in range(Nx) if y[k]==1]
            Gy = Gx.copy()

            # remove nodes not in Vy
            for Vi in Gy.nodes():
                if not Vi in Vy: Gy.remove_node(Vi)

            if len(nx.descendants(Gy, V[0])) + 1 == len(Vy):
                Ar = nx.maximum_spanning_arborescence(Gy, attr="weight")
                # assert(all(Vi in Ar.nodes() for Vi in Vy))
                # assert(len(Vy) == len(Ar.nodes()))
                logpy = sum([Ar.get_edge_data(Vi,Vj)["weight"] for Vi,Vj in Ar.edges()])
                px = max(px, exp(logpy))
        
        t += 1
    
    return px


def CALC_MAX_PROB(x, M0):
    key = TO_STRING(x[:M0])

    if key in CSET:
        CSET["hit"] += 1; return CSET[key]
        
    if CALC_LOSS(x) > L_THRES:
        CSET[key] = 1; return 1
        
    if CALC_LOSS(x[:M0]+[1]*(M-M0)) <= L_THRES:
        CSET[key] = 0; return 0

    px = VERTEX_SEARCH(x, M0)
    
    return px


# importance sampling for mu.
def EST_MU(Ns):
    Mu = [0] * Ns
    Loss = [0] * Ns
    Elap = [0] * Ns
    
    for ns in range(Ns):
        t1 = time.time()
        x, w = SAMPLE_PPLUS(ns)
        t2 = time.time()
        Mu[ns] = w
        Loss[ns] = CALC_LOSS(x)
        Elap[ns] = t2 - t1
        if DEBUG: PRINT_SIM_INFO(ns, x, w)

    return Mu, Loss, Elap


# sampling x according to the change of measure P^+
def SAMPLE_PPLUS(ns):
    x = [0] * M
    w = 1
    
    for i in range(M):
        mu_h = [0,0]

        for k in range(2):
            y = x[:i] + [k] + [0]*(M-i-1)
            mu_h[k] = CALC_MAX_PROB(y,i+1)
            
        pp = p[i] * mu_h[1] / (p[i] * mu_h[1] + (1 - p[i]) * mu_h[0])

        if rnd() < pp:
            x[i] = 1
            w *= p[i] / pp

            if mu_h[1] == 1: # no longer need to use change of measure
                for j in range(i+1,M):
                    x[j] = 1 if rnd() < p[j] else 0
                break
        else:
            x[i] = 0
            w *= (1 - p[i]) / (1 - pp)

        # if DEBUG: PRINT_SIM_INFO_DETAILED(ns, i, mu_h, pp, x)
        
    return x, w


''' brute force the value of mu_i(x_1,...,x_M0).'''
def BRUTE_FORCE(x, M0):
    if M0 == M: return (1, 1) if CALC_LOSS(x) > L_THRES else (0, 0)
    
    format = "{0:0%db}"%(M-M0)
    mu_gt = 0.0                 # optimal change of measure probability
    mu_mp = 0.0                 # minset-maxprob probability
    y = copy(x)
    
    for count in range(2**(M-M0)):
        # if count > 0 and count % 10**4 == 0:
        #     print round(100 * count / 2**(M-M0), 2)
        
        bin = format.format(count)
        
        for j in range(M-M0):
            y[M0+j] = 1 if bin[j] == "1" else 0
        
        if CALC_LOSS(y) > L_THRES:
            mu_gt += PROD([p[j] if y[j] == 1 else 1 - p[j] for j in range(M0,M)])
            mu_mp = max(mu_mp, PROD([p[j] for j in range(M0,M) if y[j] == 1]))
    
    return mu_gt, mu_mp


def NAIVE_MONTE_CARLO(Ns):
    Mu = [0] * Ns
    
    for ns in range(Ns):
        if ns % 10**4 == 0: print ns
        x = [1 if rnd() <= p[m] else 0 for m in range(M)]
        loss = CALC_LOSS(x)
        Mu[ns] = 1 if loss > L_THRES else 0

    return Mu


##################################################

# V = ['s', 'ws1', 'ws2', 'pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6', 'pc7', 'pc8', 'pc9', 'as1', 'as2', 'db1']
# L = [0, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 10, 10, 25]
# E = [('s', 'ws1'), ('s', 'ws2'), ('ws1', 'ws2'), ('ws2', 'ws1'), ('ws1', 'pc1'), ('ws1', 'pc2'), ('ws1', 'pc3'), ('ws1', 'pc4'), ('ws1', 'pc5'), ('ws1', 'pc6'), ('ws1', 'pc7'), ('ws1', 'pc8'), ('ws1', 'pc9'), ('ws1', 'as1'), ('ws1', 'as2'), ('ws2', 'pc1'), ('ws2', 'pc2'), ('ws2', 'pc3'), ('ws2', 'pc4'), ('ws2', 'pc5'), ('ws2', 'pc6'), ('ws2', 'pc7'), ('ws2', 'pc8'), ('ws2', 'pc9'), ('ws2', 'as1'), ('ws2', 'as2'), ('pc1', 'pc2'), ('pc1', 'pc3'), ('pc1', 'pc4'), ('pc1', 'pc5'), ('pc1', 'pc6'), ('pc1', 'pc7'), ('pc1', 'pc8'), ('pc1', 'pc9'), ('pc1', 'as1'), ('pc1', 'as2'), ('pc2', 'pc1'), ('pc2', 'pc3'), ('pc2', 'pc4'), ('pc2', 'pc5'), ('pc2', 'pc6'), ('pc2', 'pc7'), ('pc2', 'pc8'), ('pc2', 'pc9'), ('pc2', 'as1'), ('pc2', 'as2'), ('pc3', 'pc1'), ('pc3', 'pc2'), ('pc3', 'pc4'), ('pc3', 'pc5'), ('pc3', 'pc6'), ('pc3', 'pc7'), ('pc3', 'pc8'), ('pc3', 'pc9'), ('pc3', 'as1'), ('pc3', 'as2'), ('pc4', 'pc1'), ('pc4', 'pc2'), ('pc4', 'pc3'), ('pc4', 'pc5'), ('pc4', 'pc6'), ('pc4', 'pc7'), ('pc4', 'pc8'), ('pc4', 'pc9'), ('pc4', 'as1'), ('pc4', 'as2'), ('pc5', 'pc1'), ('pc5', 'pc2'), ('pc5', 'pc3'), ('pc5', 'pc4'), ('pc5', 'pc6'), ('pc5', 'pc7'), ('pc5', 'pc8'), ('pc5', 'pc9'), ('pc5', 'as1'), ('pc5', 'as2'), ('pc6', 'pc1'), ('pc6', 'pc2'), ('pc6', 'pc3'), ('pc6', 'pc4'), ('pc6', 'pc5'), ('pc6', 'pc7'), ('pc6', 'pc8'), ('pc6', 'pc9'), ('pc6', 'as1'), ('pc6', 'as2'), ('pc7', 'pc1'), ('pc7', 'pc2'), ('pc7', 'pc3'), ('pc7', 'pc4'), ('pc7', 'pc5'), ('pc7', 'pc6'), ('pc7', 'pc8'), ('pc7', 'pc9'), ('pc7', 'as1'), ('pc7', 'as2'), ('pc8', 'pc1'), ('pc8', 'pc2'), ('pc8', 'pc3'), ('pc8', 'pc4'), ('pc8', 'pc5'), ('pc8', 'pc6'), ('pc8', 'pc7'), ('pc8', 'pc9'), ('pc8', 'as1'), ('pc8', 'as2'), ('pc9', 'pc1'), ('pc9', 'pc2'), ('pc9', 'pc3'), ('pc9', 'pc4'), ('pc9', 'pc5'), ('pc9', 'pc6'), ('pc9', 'pc7'), ('pc9', 'pc8'), ('pc9', 'as1'), ('pc9', 'as2'), ('as1', 'pc1'), ('as1', 'pc2'), ('as1', 'pc3'), ('as1', 'pc4'), ('as1', 'pc5'), ('as1', 'pc6'), ('as1', 'pc7'), ('as1', 'pc8'), ('as1', 'pc9'), ('as1', 'as2'), ('as2', 'pc1'), ('as2', 'pc2'), ('as2', 'pc3'), ('as2', 'pc4'), ('as2', 'pc5'), ('as2', 'pc6'), ('as2', 'pc7'), ('as2', 'pc8'), ('as2', 'pc9'), ('as2', 'as1'), ('as1', 'db1'), ('as2', 'db1')]
# p = [0.1, 0.1, 0.3, 0.3, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.001, 0.001]

V = ['s', 'ws1', 'ws2', 'pc1', 'pc2', 'as1', 'as2', 'db1']
L = [0, 1, 1, 5, 5, 10, 10, 25]
E = [('s', 'ws1'), ('s', 'ws2'), ('ws1', 'ws2'), ('ws2', 'ws1'), ('ws1', 'pc1'), ('ws1', 'pc2'), ('ws1', 'as1'), ('ws1', 'as2'), ('ws2', 'pc1'), ('ws2', 'pc2'), ('ws2', 'as1'), ('ws2', 'as2'), ('pc1', 'pc2'), ('pc1', 'as1'), ('pc1', 'as2'), ('pc2', 'pc1'), ('pc2', 'as1'), ('pc2', 'as2'), ('as1', 'pc1'), ('as1', 'pc2'), ('as1', 'as2'), ('as2', 'pc1'), ('as2', 'pc2'), ('as2', 'as1'), ('as1', 'db1'), ('as2', 'db1')]
p = [0.1, 0.1, 0.3, 0.3, 0.01, 0.01, 0.001, 0.001, 0.01, 0.01, 0.001, 0.001, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.001, 0.001]

# V = ['s', 'ws1', 'ws2', 'pc1', 'as1', 'as2', 'db1']
# L = [0, 1, 1, 5, 10, 10, 25]
# E = [('s', 'ws1'), ('s', 'ws2'), ('ws1', 'ws2'), ('ws2', 'ws1'), ('ws1', 'pc1'), ('ws1', 'as1'), ('ws1', 'as2'), ('ws2', 'pc1'), ('ws2', 'as1'), ('ws2', 'as2'), ('pc1', 'as1'), ('pc1', 'as2'), ('as1', 'pc1'), ('as1', 'as2'), ('as2', 'pc1'), ('as2', 'as1'), ('as1', 'db1'), ('as2', 'db1')]
# p = [0.1, 0.1, 0.3, 0.3, 0.01, 0.001, 0.001, 0.01, 0.001, 0.001, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.001, 0.001]

M = len(E)
N = len(V)

logp = [log(pi) for pi in p]
ratio = 0.9 #float(sys.argv[1])
L_MAX = sum(L)
L_THRES = ratio * L_MAX
Ns = 1000
DEBUG = True
SIM_INFO = "M=%d N=%d Ns=%s L_MAX=%.1f ratio=%.2f"%(M, N, Ns, L_MAX, ratio)
CSET = {"info":SIM_INFO, "hit":0}

assert(V[0] == "s")
assert(len(L) == N)
assert(len(p) == M)

print SIM_INFO
print "L=", L
print "p=", p

##################################################

# EST_MU(1)
# VERTEX_SEARCH([0]*M,1)
# VERTEX_SEARCH([1]+[0]*(M-1),0)
# exit()

# filename = "campus_N_%d_ratio_%d.txt"%(N, ratio*100)

# while True:
#     Mu, Loss, Elap = EST_MU(10)
#     openfile = open(filename, "a")
#     openfile.write(",".join(map(str,Mu)) + "\n")
#     openfile.close()
    
# exit()

Mu, Loss, Elap = EST_MU(Ns)
mu_mean = np.mean(Mu)
print SIM_INFO
print "mu_mean", mu_mean

for k in [1, 50, 100, 500, 1000]:
    Mu_k = AVERAGE_EVERY(Mu, k)
    print k, CALC_RELATIVE_ERROR(Mu_k)

exit()
