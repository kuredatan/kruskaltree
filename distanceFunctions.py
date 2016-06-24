from graph import Graph

#Distance function
#Distance functions take as argument the super tree (tree of all trees), and returns a float resulting from the comparison of the two nodes of the super tree
#node are integers (identifiers for TaxoTrees. Conversion to TaxoTree is allowed by @hashArray in the graph)

def printDistance():
    print "LCAdifference"
    return ["LCAdifference"]

def applyFctD(fct,superTree,node1,node2):
    if fct == "LCAdifference":
        return LCAdiff(superTree,node1,node2)

#Node indexed with 0 is arbitrarily considered the root of the SuperTree
#Starting from node and going up to root
def getPath(superTree,node):
    #Resulting path: list of integers from root to node
    path = [node]
    currNode = node
    while not (currNode == 0):
        nextEdge = [ i for (i,j,k) in superTree.edges if j == node ]
        if not (len(nextEdge) == 1):
            print "\n/!\ ERROR: Not a tree."
            raise ValueError
        path = [nextEdge[0]] + path
        if currNode == nextEdge[0]:
            print "\n/!\ ERROR: Infinite loop."
            raise ValueError
        else:
            currNode = nextEdge[0]
    return path

#Returns common path
def getCommon(path1,path2):
    commonPath = []
    #path1 must be at least of length 1 (for root)
    if not path1 or not path2:
        print "\n/!\ ERROR: Empty path."
        raise ValueError
    currNode1 = path1.pop()
    currNode2 = path2.pop()
    while path1 and path2 and currNode1 == currNode2:
        commonPath.append(currNode1)
        currNode1 = path1.pop()
        currNode2 = path2.pop()
    #If both @path1 and @path2 are empty here, it means
    #@commonPath = @path1 = @path2
    if path1:
        path1.append(currNode1)
    if path2:
        path2.append(currNode2)
    return commonPath

#returns |total path1| + |total path2| - 2*|common path|
def LCAdiff(superTree,node1,node2):
    path1 = getPath(superTree,node1)
    path2 = getPath(superTree,node2)
    commonPath = getCommon(path1,path2)
    n1 = len(path1)
    n2 = len(path2)
    return n1 + n2
