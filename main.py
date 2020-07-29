from __future__ import division
import networkx as nx
from copy import copy
import random
from random import random as rnd, randint as rndi
from math import log, log10, exp, sqrt
import matplotlib.pyplot as plt
import numpy as np
import time

from util import *

# Configure save dialog to use current directory
import matplotlib as mpl
mpl.rcParams["savefig.directory"] = ""


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

    # Vx2 contains all vertices that may be reachable from s
    Vx2 = list(nx.descendants(Gx, V[0]))
    assert(all(Vi in Vx2 for Vi in Vx1))
    
    # Remove nodes in Gx that are not in Vx2
    for Vi in Gx.nodes():
        if Vi != V[0] and not Vi in Vx2: Gx.remove_node(Vi)
    
    # Initial setup for vertex search
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

            # Remove nodes not in Vy
            for Vi in Gy.nodes():
                if not Vi in Vy: Gy.remove_node(Vi)

            if len(nx.descendants(Gy, V[0])) + 1 == len(Vy):
                Ar = nx.maximum_spanning_arborescence(Gy, attr="weight")
                logpy = sum([Ar.get_edge_data(Vi,Vj)["weight"] for Vi,Vj in Ar.edges()])
                px = max(px, exp(logpy))
        
        t += 1
    
    return px


def CALC_MINSET_MAXPROB(x, M0):
    key = TO_STRING(x[:M0])

    if key in CSET:
        CSET["hit"] += 1; return CSET[key]
        
    if CALC_LOSS(x) > L_THRES:
        CSET[key] = 1; return 1
        
    if CALC_LOSS(x[:M0]+[1]*(M-M0)) <= L_THRES:
        CSET[key] = 0; return 0

    px = VERTEX_SEARCH(x, M0)
    CSET[key] = px
    
    return px


# Importance sampling for mu.
def EST_MU(Ns):
    Mu = [0] * Ns
    Loss = [0] * Ns
    Elap = [0] * Ns
    
    for ns in range(Ns):
        start_time = time.time()
        x, w = SAMPLE_PPLUS(ns)
        end_time = time.time()
        Mu[ns] = w
        Loss[ns] = CALC_LOSS(x)
        Elap[ns] = end_time - start_time
        
        if DEBUG: PRINT_SIM_INFO(ns, x, w)
        
    return Mu, Loss, Elap


# Sample x according to the change of measure P^+
def SAMPLE_PPLUS(ns):
    x = [0] * M
    w = 1
    
    for i in range(M):
        mu_h = [0, 0]

        for k in range(2):
            y = x[:i] + [k] + [0] * (M-i-1)
            mu_h[k] = CALC_MINSET_MAXPROB(y, i+1)
            
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


''' Brute force the value of mu_i(x_1,...,x_M0).'''
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


def COMPLETE_STATE_ENUMERATION():
    format = "{0:0%db}"%M
    Loss = {}
    
    for count in range(2**M):
        bin = format.format(count)

        y = [1 if bin[j] == "1" else 0 for j in range(M)]

        ly = CALC_LOSS(y)
        py = PROD([p[j] if y[j] == 1 else 1 - p[j] for j in range(M)])

        if not ly in Loss: Loss[ly] = py
        else: Loss[ly] += py

    return Loss

    
def NAIVE_MONTE_CARLO(Ns):
    Mu = [0] * Ns
    
    for ns in range(Ns):
        if ns % 10**4 == 0: print ns
        x = [1 if rnd() <= p[m] else 0 for m in range(M)]
        loss = CALC_LOSS(x)
        Mu[ns] = 1 if loss > L_THRES else 0

    return Mu


##################################################

from campus1 import *

M = len(E)
N = len(V)

logp = [log(pi) for pi in p]
ratio = 0.99 #float(sys.argv[1])
L_MAX = sum(L)
L_THRES = ratio * L_MAX
Ns = 10000
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

filename = "campus_N_%d_ratio_%d.txt"%(N, ratio*100)

while True:
    Mu, Loss, Elap = EST_MU(10)
    openfile = open(filename, "a")
    openfile.write(",".join(map(str,Mu)) + "\n")
    openfile.close()
    
exit()
