#Complete graph of all trees in adjacency lists
#Stores trees as integers (hashArray[i] is the unique tree associated to integer i) and edges as 3-ples of integers (first node, second node,distance)
#For superTree (tree of all trees obtained with Kruskal's algorithm) root will be chosen arbitrarily

from taxoTree import TaxoTree
from comparisonFunctions import applyFctC

class Graph(object):
    #@n is the total number of trees/patients
    def __init__(self,n=0,vertices=None,edges=None,hashArray=None):
        self.vertices = vertices or []
        self.edges = edges or []
        self.hashArray = hashArray or [None]*n
        #
        #
    #f is a comparison function: takes as arguments two sample names groups and the whole TaxoTree, gives a non-negative weight
    #Gives the complete graph of all treeswith weighted edges according to the distance function f provided
    def constructComplete(self,samplesGroupList,tree,f):
        n = len(samplesGroupList)
        hashArray = []
        #Initializing the hashArray
        for samplesGroup in samplesGroupList:
            hashArray.append(samplesGroup)
        self.hashArray = hashArray
        self.vertices = [i for i in range(n)]
        #The orientation of the graph depends on f
        self.edges = [(i,j, applyFctC(f,self.hashArray[i],self.hashArray[j])) for i in range(n) for j in range(n) if not (j == i)]
        return self
    #
    #
    #When Kruskal's algorithm is applied to the complete graph, one creates a new graph, at first empty, to which edges are added one by one
    #completeHashArray is the hashArray of the complete graph
    def addVertex(self,completeHashArray,edge):
        i,j,w = edge
        n = len(completeHashArray)
        #Allows to get the trees in growing order of the hash integers
        #"not hashArray[k]" means "not (hashArray[k] == None)" ie. hashArray[k] contains a tree that belongs to the graph
        newHash = []
        for k in range(n):
            if (k == i) or (k == j) or (k < n and not self.hashArray[k]):
                newHash.append(completeHashArray[k])
            else:
                newHash.append(None)
        self.hashArray = newHash
        self.vertices = [i for i in range(n) if not self.hashArray[i]]
        self.edges.append(edge)
        return self
