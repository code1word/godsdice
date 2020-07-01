# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 17:15:34 2020
@author: Thoma
"""

# %% codecell
# imports

import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
from numpy import kron
from qiskit.visualization import plot_histogram

# %% codecell
# setup functions

#Adjacency standard of Aij is a line from i to j

def indegrees(M,node):
    #This is the sum of the colume at node X, counting from 0
    return sum(M[:,node])

def outdegrees(M,node):
    #Counts the entries in the node row
    return sum(M[node])

def aVal(M,x):
    indegree = indegrees(M,x)
    outdegree = outdegrees(M,x)
    if indegree == 0 and outdegree == 0:
        return 1
    else:
        return indegree / (indegree + outdegree)

def constructC(M):

    # Force input matrix to already have dimensions as a power of 2
    if not np.isclose(np.log2(M.shape[0]), int(np.log2(M.shape[0]))):
        print("Please pad your adjacency matrix so its dimensions are a power of 2")

    output = np.zeros([2 * x for x in M.shape])

    for x in range(M.shape[0]):

        xVec = np.zeros(M.shape[0])
        xVec[x] = 1
        outerX = np.diag(xVec)

        a = aVal(M,x)
        cellElement = np.array([[np.sqrt(1/(a + 1)),np.sqrt(a/(1 + a))],
                                [np.sqrt(a/(1 + a)),-1 * np.sqrt(1/(a + 1))]
                                ])
        output += kron(cellElement,outerX)

    return output

def constructS(M):

    # Force input matrix to already have dimensions as a power of 2
    if not np.isclose(np.log2(M.shape[0]), int(np.log2(M.shape[0]))):
        print("Please pad your adjacency matrix so its dimensions are a power of 2")

    up_proj = np.diag([1,0])
    down_proj = np.diag([0,1])

    P,D,Q = la.svd(M)
    U = P @ scila.expm(1j*np.diag(D)) @ Q
    if_up = np.kron(up_proj, np.identity(M.shape[0]))
    if_down = np.kron(down_proj, U)
    S = if_up+if_down

    return S

# adjacency matrix for binary trees
def binaryTree(levels):
    n_nodes = 2**(levels+1)
    # final node recieving a connection
    n_final = 2**levels - 1

    adj_mat = np.zeros([n_nodes, n_nodes])
    for j in range(n_final):
        adj_mat[2*j+1, j] = adj_mat[2*j+2, j] = 1

    return adj_mat

# adjacency matrix for scale-free network
def scaleFree(nodes):
    adj_mat = np.zeros([nodes,nodes])

    adj_mat[0,0]=1

    for node in range(1,nodes):
        probs = np.array([sum(adj_mat[i,:])+sum(adj_mat[:,i]) for i in range(node)])
        if not sum(probs) == 0:
            probs = probs/sum(probs)

        for join in range(node):
            if np.random.random() < probs[join]:
                adj_mat[node, join] = 1
    return adj_mat

# %% codecell
# 7-node graph from paper
A = np.array([  [0,1,0,0,1,1,1,0],
                [0,0,0,0,0,0,0,0],
                [1,1,0,0,0,0,1,0],
                [0,0,1,0,1,1,0,0],
                [0,0,0,0,0,0,1,0],
                [0,0,1,0,0,0,0,0],
                [0,0,0,0,1,0,0,0],
                [0,0,0,0,0,0,0,0]])

# 5-level binary tree
# A = binaryTree(5)

# 6-qubit scale-free network
# A = scaleFree(32)

# initialize statevector, operators
C = constructC(A)
S = constructS(A)
SC = S @ C

initialState = np.zeros(8)
active_nodes = 7

for i in range(active_nodes):
    initialState[i] = 1/np.sqrt(active_nodes)

psio = kron(np.array([1,1]),initialState)/np.sqrt(2)

# %% codecell
# run the circuit - get SV after one step, then take one more step, etc.
finalStates = [SC@psio]
for i in range(500):
    finalStates.append(SC @ finalStates[i])

# find probabilities (instantaneous rank) after each step
probabilityList = []
for item in range(len(finalStates)):
    probabilities = [np.real(p * np.conj(p)) for p in finalStates[item]]
    probabilityList.append(probabilities)

# combine instantaneous ranks into final average rank
# creates as many dictionary slots as we need
counts = {}
for i in range(2**3):
    counts[bin(i)[2:].zfill(3)] = 0

for item in probabilityList[50:]:
    for i in range(len(item)):
        # figures out the binary label that qiskit would return
        binI = bin(i)[2:].zfill(4)[1:]
        counts[binI] += item[i]

plot_histogram(counts)
