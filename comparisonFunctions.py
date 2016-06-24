from __future__ import division

from taxoTree import TaxoTree
from misc import takeNodesInTree

#Comparison with total Ratio, pattern Ratio, tree edit distance (e.g. Zhang-Shasha's algorithm), etc.
#Comparison functions take as argument two groups of samples/directly the tree associated to the groups, computes if needed the (forest of) trees associated to these groups, and returns a float resulting from the comparison of the trees

def printComparison():
    print "f"
    return ["f"]

#Will return a list of TaxoTrees associated to @sampleGroup
#Complexity?
def convertSampleIntoTree(tree,sampleGroup):
    return [TaxoTree("TODO")]

def applyFctC(fct,tree,sampleGroup1,sampleGroup2):
    if fct == "f":
        return f(convertSampleIntoTree(tree,sampleGroup1),convertSampleIntoTree(tree,sampleGroup2))

#(Dumb) test function
def f(tree1,tree2):
    return 1
