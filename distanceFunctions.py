from __future__ import division

#Distance function
#Distance functions take as argument the super tree (tree of all trees), and returns a float resulting from the comparison of the two nodes of the super tree

def printDistance():
    print "LCAdifference"
    return ["LCAdifference"]

def applyFctD(fct,superTree,node1,node2):
    if fct == "LCAdifference":
        return LCAdiff(superTree,node1,node2)

#convert the graph into a tree
    
#the node indexed by 0 is considered root of the super tree

#returns |path1| + |path2| - 2*|common path|
def LCAdiff(superTree,node1,node2):
    return 1
