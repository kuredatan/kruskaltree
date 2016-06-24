import numpy as np
import re

from writeOnFiles import writeFile
from taxoTree import TaxoTree,printTree
from misc import mem,isInDatabase,partitionSampleByMetadatumValue,setOperations

from graph import Graph
from kruskal import kruskal
from comparisonFunctions import printComparison,f
from distanceFunctions import printDistance,applyFctD

#@dataArray = [samplesInfoList,infoList,samplesOccList,speciesList,paths,n,nodesList,taxoTree,sampleIDList]

integer = re.compile("[0-9]+")

#Parsing functions
def parseList(string):
    if not (len(string.split(",")) == 1):
        print "\n/!\ ERROR: Do not use ',' as a separator: rather use ';'."
        raise ValueError
    elif not (len(string.split(":")) == 1):
        print "\n/!\ ERROR: Do not use ':' as a separator: rather use ';'."
        raise ValueError
    return string.split(";")

def parseListNode(string):
    if not (len(string.split(":")) == 1):
        print "\n/!\ ERROR: Do not use ':' as a separator: rather use ';'."
        raise ValueError
    ls = string.split(";")
    res = []
    for node in ls:
        nodeSplit = node.split(",")
        if not (len(nodeSplit) == 2):
            print "\n/!\ ERROR: Please use ',' as a separator for name,rank of a bacteria."
            raise ValueError
        nodeSplitName = nodeSplit[0].split("(")
        if not (len(nodeSplitName) == 2):
            print "\n/!\ ERROR: Please use the syntax '([name],[rank])' for each bacteria."
            raise ValueError
        nodeSplitRank = nodeSplit[-1].split(")")
        if not (len(nodeSplitRank) == 2):
            print "\n/!\ ERROR: Please use the syntax '([name],[rank])' for each bacteria."
            raise ValueError
        name,rank = nodeSplitName[-1],nodeSplitRank[0]
        res.append((name,rank))
    return res

def parseIntList(string):
    if not (len(string.split(",")) == 1):
        print "\n/!\ ERROR: Do not use ',' as a separator: rather use ';'."
        raise ValueError
    elif not (len(string.split(":")) == 1):
        print "\n/!\ ERROR: Do not use ':' as a separator: rather use ';'."
        raise ValueError
    l = string.split(";")
    resultList = []
    for s in l:
        if integer.match(s):
            resultList.append(int(s))
        elif s == "+inf" or s == "-inf":
            resultList.append(s)
        else:
            print "\n/!\ ERROR: Here you can only use integers or '+inf' or '-inf'."
            raise ValueError
    return resultList

#___________________________________________________________________________

#Macros for formatting
#Printing pretty lists of nodes
def listNodes(nodeList):
    string = ""
    for l in nodeList[:-1]:
        string += str(l) + ", "
    string += str(nodeList[-1])
    return string

#@stringNode is assumed to be a (name,rank) pair, with name and rank being strings
#@sanitizeNode allows it to be printed "(name,rank)" and not "('name','rank')"
def sanitizeNode(stringNode):
    return "(" + stringNode[0] + "," + stringNode[1] + ")"

#Printing pretty lists of metadata with their default values
def listSampleInvolved(metadataList,interval1List,interval2List,sampleNameList):
    string = ""
    if not metadataList and not interval1List and not interval2List and not sampleNameList:
      print "\n/!\ ERROR: You have selected no sample."
      raise ValueError
    #If samples were selected one by one
    elif sampleNameList:
        string += "\ndepending on the group of samples: "
        for sl in sampleNameList[:-1]:
            string += str(sl) + ", "
        string += str(sampleNameList[-1])
    #If samples were selected according to metadata values (len(metadataList) = len(interval1List) = len(interval2List))
    if metadataList:
        string += "\nselected on metadata (for each line): "
        n = len(metadataList)
        for i in range(n-1):
            if (interval1List[i] == interval2List[i]):
                string += metadataList[i] + " (value equal to " + str(interval1List[i]) + "), "
            else:
                string += metadataList[i] + " (value between " + str(interval1List[i]) + " and " + str(interval2List[i]) + "), "
        if (interval1List[-1] == interval2List[-1]):
            string += metadataList[-1] + " (value equal to " + str(interval1List[-1]) + ")"
        else:
            string += metadataList[-1] + " (value between " + str(interval1List[-1]) + " and " + str(interval2List[-1]) + ")"
    return string

#Selecting samples in two ways: either choose each of them one by one, or selecting according to default values of certain metadatum
def createSampleNameList(dataArray):
    metadataList = []
    interval1List = []
    interval2List = []
    sampleIDList = dataArray[8]
    i = raw_input("/!\ How many different lists of samples do you want?\n")
    if not integer.match(i):
        print "\n/!\ ERROR: You need to enter a integer here!"
        raise ValueError
    numberList = int(i)
    sampleNameList = []
    if (numberList < 1):
        print "\n/!\ ERROR: Empty set of lists of samples!"
        raise ValueError
    while numberList:
        answer = raw_input("Do you want to select samples one by one, or to select samples matching requirements on metadata? one/matching \n")
        if (answer == "one"):
            if (len(sampleIDList) < 2):
                print "\n/!\ ERROR: List of samples is empty or only of length one!..."
                raise ValueError
            print sampleIDList
            sampleNameList11 = parseList(raw_input("Input the list of samples using the ID printed above. [e.g. " + sampleIDList[0] + ";"+ sampleIDList[1] + " ]\n"))
        elif (answer == "matching"):
            print dataArray[1]
            metadataList = parseList(raw_input("Input the list of metadata you want to consider among those written above. [ e.g. " + dataArray[1][0] + ";" + dataArray[1][-1] + " ]\n"))
            isInDatabase(metadataList,dataArray[1])
            interval1List = parseIntList(raw_input("Input the list of lower interval bounds corresponding to metadatum/metadata above. [ Please refer to README for more details. e.g. 1;2 ]\n"))
            if not (len(interval1List) == len(metadataList)):
                print "\n/!\ ERROR: You need to enter the same number of lower bounds than of metadata!"
                raise ValueError
            interval2List = parseIntList(raw_input("Input the list of upper interval bounds corresponding to metadatum/metadata above. [ Please refer to README for more details. e.g. 3;2 ]\n"))
            if not (len(interval2List) == len(metadataList)):
                print "\n/!\ ERROR: You need to enter the same number of upper bounds than of metadata!"
                raise ValueError
            sampleNameList11 = computeSamplesInGroup(dataArray[0],dataArray[1],metadataList,interval1List,interval2List)[0]
        else:
            print "\n/!\ ERROR: You need to answer either 'one' or 'matching' and not: \"",answer,"\"."
            raise ValueError
        isInDatabase(sampleNameList11,sampleIDList)
        sampleNameList.append(sampleNameList11)
        numberList -= 1
    return sampleNameList,metadataList,interval1List,interval2List

#____________________________________________________________________________

#Actions
def runAct(dataArray):
    print "Choosing the list of samples."
    #or use partition by metadatum values
    sampleNameList,metadataList,interval1List,interval2List = createSampleNameList(dataArray)
    n = len(sampleNameList)
    print "\nAVAILABLE COMPARISON FUNCTION(S):"
    fctF = printComparison()
    f = raw_input("\nChoose your comparison function above those printed above.\n")
    isInDatabase([f],fctF)
    completeGraph = Graph(n).constructComplete(sampleNameList,dataArray[7],f)
    superTree,w = kruskal(completeGraph)
    #Constructing distance matrix
    matrix = np.zeros((n,n))
    print "\nAVAILABLE DISTANCE FUNCTION(S):"
    fctD = printDistance()
    d = raw_input("\nChoose your distance function above those printed above.\n")
    isInDatabase([d],fctD)
    valueArray = []
    print "\nSUPERTREE of weight:",w
    print superTree.vertices
    print superTree.edges
    for i in range(n):
        for j in range(i,n):
            #matrix is symmetric (distance)
            s = applyFctD(d,superTree,i,j)
            matrix[i][j] = s
            matrix[j][i] = s
            valueArray.append(s)
    valueArray = sorted(valueArray)
    valueNumber = n*n/2
    quartile3 = valueNumber*3/4
    valueQuartile = valueArray[quartile3]
    mostDifferent = []
    #Distance is symmetric
    for i in range(n):
        for j in range(i+1,n):
            if matrix[i][j] >= valueQuartile:
                mostDifferent.append((sampleNameList[i],sampleNameList[j]))
    print "\nRESULTING MATRIX:"
    print matrix
    print "\n---\nMost different samples groups from:\n"
    for sampleGroup in sampleNameList:
        print sampleGroup
    print "\nare:\n"
    print mostDifferent
    print "\n--- END OF DISPLAY\n"

#____________________________________________________________________________

def printTreeAct(dataArray):
    answer = raw_input("Do you want to print sample hit lists? Y/N\n")
    if not ((answer == "Y") or (answer == "N")):
        print "\n/!\ ERROR: You need to answer 'Y' or 'N'."
        raise ValueError
    printTree(dataArray[7],(answer == "Y"))
