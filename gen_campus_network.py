from __future__ import division

    
def add_edge(V1, V2, level):
    global E, p
    E += [(V1, V2)]
    p += [THREAT[level]]

    
def add_node(Vi, value):
    global V, L
    V += [Vi]
    L += [LOSS[value]]


S = ["s"]                       # initial point of intrusion
WS = ["ws1", "ws2"]             # web servers in DMZ
U = ["pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9"] # personal computers in Internal Network
AS = ["as1", "as2"]              # application servers in Internal Network
DB = ["db1"]                     # database server in Office Network

THREAT = {"certain":1.0, "highly-likely":0.5, "probable":0.3, "likely":0.1,
          "moderate":0.03, "unlikely":0.01, "remote":1e-3, "highly-unlikely":1e-4,
          "impossible":0.0}

LOSS = {"none":0, "small":1, "moderate":5, "high":10, "critical":25}


V, L, E, p = [], [], [], []


for Vi in S + WS + U + AS + DB:
    if Vi in S: add_node(Vi, "none")
    if Vi in WS: add_node(Vi, "small")
    if Vi in U: add_node(Vi, "moderate")
    if Vi in AS: add_node(Vi, "high")
    if Vi in DB: add_node(Vi, "critical")


for V1 in S:
    for V2 in WS:
        add_edge(V1, V2, "likely")

for V1 in WS:
    for V2 in WS:
        if V1 != V2: add_edge(V1, V2, "probable")
            
for V1 in WS:
    for V2 in U:
        add_edge(V1, V2, "unlikely")

    for V2 in AS:
        add_edge(V1, V2, "remote")

for V1 in U + AS:
    for V2 in U + AS:
        if V1 != V2: add_edge(V1, V2, "moderate")

for V1 in AS:
    for V2 in DB:
        add_edge(V1, V2, "remote")

        
for i in range(len(V)):
    print i, L[i], V[i]

    
for i in range(len(E)):
    print i, p[i], E[i]

    
print "V =", V
print "L =", L
print "E =", E
print "p =", p
print len(V), len(E)
