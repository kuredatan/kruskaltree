from __future__ import division
import numpy as np
import re

from writeOnFiles import writeFile

integer = re.compile("[0-9]+")

#Returns a list of (name,rank,sampleHitList) tuples associated with
#nodes of the taxonomic tree reduced to the sample sampleName
#and the number of assignments in this tree
#NB: We need to keep the whole sampleHitList (and not only the
#number associated with sampleName) in order to apply set operations (intersection, ...)
def inSample(element,sampleNameList):
    #element[1] (number associated to a sample) must be non-zero
    if not element[1]:
        print "\n/!\ ERROR: [BUG] [misc/inSample] Null element."
        raise ValueError
    #If list is empty
    if not sampleNameList:
        return False
    for name in sampleNameList:
        if (element[0] == name):
            return True
    return False

#@sampleName is a list of names of samples
#takeNodesInTree(tree,sampleName) should return the taxonomic tree reduced to the nodes assigned in sample sampleName = the list of assigned node + their sampleHitList, and the number of assignments in the reduced tree, and the number of nodes in the reduced tree
def takeNodesInTree(tree,sampleNameList):
    #@sampleHitList (see TaxoTree) is a list of (name of sample, number) pairs attached to a node of a TaxoTree
    sample = []
    numberTotalAssignments = 0
    numberNodes = 0
    queue = [ tree ]
    while queue:
        node = queue.pop()
        isInSample = []
        for x in node.sampleHitList:
            if inSample(x,sampleNameList):
                isInSample.append(x)
        #if node is in sample, ie isInSample is not empty
        if isInSample:
            sample.append((node.name,node.rank,node.sampleHitList))
            numberTotalAssignments += isInSample[0][1]
            numberNodes += 1
        queue += node.children
    return sample,numberTotalAssignments,numberNodes

#Returns boolean and sampleHitList if true
def memAndSampleHitList(x,nodeList):
    sampleHitList = []
    nodeListCopy = []
    for nd in nodeList:
        nodeListCopy.append(nd)
    #While @nodeList is not empty and @sampleHitList is empty
    while nodeListCopy:
        node = nodeListCopy.pop()
        if (x[0] == node[0] and x[1] == node[1]):
            return True,node[2]
    return False,[]

#Gets sample IDs from the data matrix
#/!\ Some of the samples may be appear in the data matrix!
def getSampleIDList(samplesList):
    sampleIDList = []
    for sample in samplesList:
        if not mem(sample[0],sampleIDList):
            sampleIDList.append(sample[0])
    #Sorts sample IDs in alphabetical order
    return sorted(sampleIDList,key=lambda x:x)

#For sorting lists of nodes (name,rank) by decreasing order S > G > F > O > C > P > K > R. Returns 1 if rank1 => rank2, -1 if rank1 < rank2
#You can modify the ranks and the order by modifying the default rank array
def compare(rank1,rank2,ranks=['S','G','F','O','C','P','K','R'],ranksArrayLength=8):
    for i in range(ranksArrayLength):
        if (rank2 == ranks[i]):
            return 1
        elif (rank1 == ranks[i]):
            return -1
    print "\n/!\ ERROR: Unknown rank. Did you modify the default rank array?"
    raise ValueError

def sanitize(name):
    ls = name.split(" ")
    if (len(ls) == 1):
        return ls[0]
    sName = ""
    sLs = []
    for l in ls:
        if not (l == "" or l == "(class)" or l == "\n" or l == "#"):
            sLs.append(l)
    for l in sLs[:-1]:
        sName = sName + l + " "
    sName = sName + sLs[-1]
    return sName.split("\n")[0]

#is member function
def mem(x,ls):
    n = len(ls)
    for i in range(n):
        if (x == ls[i]):
            return True
    return False

def containsSpecie(path,name,rank):
    for x in path:
        if (x[0] == name) and (x[1] == rank):
            return True
    return False

#@paths is the paths list of a TaxoTree
#@n = len(@paths)
def selectPath(paths,name,rank,n):
    i = 0
    while i < n and not containsSpecie(paths[i],name,rank):
	i += 1
    if (i == n):
        print "\n/!\ ERROR: [BUG] [misc/selectPath] No path for (%s,%s)."%(name,rank)
        raise ValueError
    else:
        path = []
        x = paths[i]
        #Memorizes only the part of the path that leads to (name,rank)
        n = len(x)
        i = 0
        while (i < n) and not (x[i][0] == name and x[i][1] == rank):
            path.append(x[i])
            i += 1
        return path

def setOperations(paths,name1,rank1,name2,rank2,allNodes=False):
    path1 = selectPath(paths,name1,rank1)
    path2 = selectPath(paths,name2,rank2)
    n = min(len(path1),len(path2))
    #if there is more than one path to the nodes, or no path
    if (n < 1):
        print "\n/!\ ERROR: [BUG] [misc/setOperations] It is not a tree."
        raise ValueError
    else:
        commonPath = []
        i = 0
        #As long as path1 and path2 are not empty
        while (i < n) and (path1[i] == path2[i]):
            commonPath.append(path1[i])
            i += 1
        return commonPath,path1[i+1:],path2[i+1:]

#Computes LCA from the list paths of a TaxoTree
def taxoLCA(paths,name1,rank1,name2,rank2,allNodes=False):
    commonPath,_,_ = setOperations(paths,name1,rank1,name2,rank2,allNodes=False)
    return commonPath[-1]

#Checks if the elements in @parselist belong to @datalist else returns an error
def isInDatabase(parseList,dataList):
    for pl in parseList:
        if not mem(pl,dataList):
            n = len(dataList)
            if not n:
                print "\n/!\ ERROR: [BUG] [actions/isInDatabase] Empty list."
            else:
                print "\n/!\ ERROR: '" + str(pl) + "' is not in the database beginning with: " + str(dataList[:min(n-1,3)]) + "."
            raise ValueError

#Given a set of samples, gives the list of disjoint groups of samples according to the value of the metadatum, and the set of values of the metadatum
#@metadatum is a list (of one element) of metadata.
def partitionSampleByMetadatumValue(metadatum,infoList,samplesInfoList):
    #One metadatum only!
    metadatum = metadatum[0]
    #computes the number of column which matches the metadatum in infoList
    i = 0
    n = len(infoList)
    while i < n and not (infoList[i] == metadatum):
        i += 1
    if (i == n):
        print "\n/!\ ERROR: metadatum",metadatum,"not found"
        raise ValueError
    #Getting the set of values of the metadatum
    #Sorting samples according to the values of the metadatum
    sampleSorted = sorted(samplesInfoList,key=lambda x: x[i])
    #List of list of samples: one sublist matches a value of the metadatum
    valueSampleMetadatum = []
    #The set of values of the metadatum
    valueSet = []
    if not len(sampleSorted):
        print "\n/!\ ERROR: You have selected no sample."
        raise ValueError
    sample = sampleSorted.pop()
    if len(sample) < i:
        print "\n/!\ ERROR: [BUG] [misc/partitionSampleByMetadatumValue] Different lengths",len(sample),"and",i,"(1)"
        raise ValueError
    #selects a sample where the value of the metadatum is known
    while not integer.match(sample[i]):
        sample = sampleSorted.pop()
        if len(sample) < i:
            print "\n/!\ ERROR: [BUG] [misc/partitionSampleByMetadatumValue] Different lengths",len(sample),"and",i,"(2)"
            raise ValueError
    #Initializing the set of values of the metadatum
    currValue = sample[i]
    valueSet.append((metadatum,int(currValue)))
    #While it remains samples in the list
    while sampleSorted:
        valueSample = []
        isEmpty = False
        #Filling the list of samples with similar values of the metadatum
        while sampleSorted and (sample[i] == currValue):
            valueSample.append(sample)
            sample = sampleSorted.pop()
            #gets the next sample where the value of the metadatum is known
            while not integer.match(sample[i]) and sampleSorted:
                sample = sampleSorted.pop()
            #If sampleSorted is empty
            if not sampleSorted:
                isEmpty = True
        #appends the newly created list to the main list
        valueSampleMetadatum.append(valueSample)
        #Initializing next loop with the new different value of the metadatum
        currValue = sample[i]
        if isEmpty:
            #Adding this value to the set
            valueSet.append((metadatum,int(currValue)))
    return valueSet,valueSampleMetadatum
