#Implementation of Kruskal using Union Find
from taxoTree import TaxoTree

from graph import Graph
import numpy as np
#UnionFind comprises fathers and ranks arrays: father[i] is the number of i's father and ranks[i] is the estimated height of the tree/class rooted at i

def initUnionFind(n):
    fathers = np.zeros(n)
    ranks = np.zeros(n)
    for i in range(n):
        fathers[i] = i
        ranks[i] = 1
    return fathers,ranks

#Find the representant of i's class:
#if i's father is the highest ancestor (ie. i's grandfather is its father), then i's father is the representant
#else we apply the algorithm to i's father
#Path compression: at the end of the search, every visited node is linked to root
def find(fathers,ranks,i):
    currNode = i
    nodeToRoot = [i]
    while not (fathers[currNode] == currNode):
        currNode = fathers[currNode]
        nodeToRoot.append(currNode)
    #Minus one because the second to last visited node is already linked to root
    n = len(nodeToRoot) - 1
    ranks[currNode] = ranks[i] - n
    for node in nodeToRoot:
        fathers[i] = currNode
    return fathers,ranks,currNode

#Implements path compression
def union(fathers,ranks,i,j):
    fathers,ranks,repi = find(fathers,ranks,i)
    fathers,ranks,repj = find(fathers,ranks,j)
    if repi == repj:
        return fathers,ranks,False
    if (ranks[i] < ranks[j]):
        fathers[repj] = repi
        ranks[repi] = ranks[repj] + 1
    else:
        fathers[repi] = repj
        ranks[repj] = ranks[repi] + 1
    return fathers,ranks,True

def compare(e1,e2):
    i,j,k = e1
    l,m,q = e2
    return (k < q)

def kruskal(graph):
    n = len(graph.vertices) + 1
    superTree = Graph(n-1,graph.vertices)
    sortedEdges = sorted(graph.edges,key = lambda e: e[2])
    fathers,ranks = initUnionFind(n)
    weight = 0
    for (i,j,w) in sortedEdges:
        fathers,ranks,boolean = union(fathers,ranks,i,j)
        if boolean:
            weight += w
            superTree = superTree.addEdge(graph.hashArray,(i,j,w))
    return (superTree,weight)


#TESTS:
def f(tree1,tree2):
    i = int(tree1.name)
    j = int(tree2.name)
    if i<j:
        return (i+10)
    else:
        return (j+1)

def kruskalTest(n=3):
    vertices = [i+1 for i in range(n)]
    hashArray = [TaxoTree("%d"%i) for i in range(n)]
    completeGraph = Graph(1).constructComplete(hashArray,f)
    print completeGraph.edges
    t,w = kruskal(completeGraph)
    print t.edges,w
