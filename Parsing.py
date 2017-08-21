#!/bin/bash
import sys
import math
import decimal
from collections import defaultdict

import Queue



def buildDicts(data):
    data = open(data, 'r')
    rhsDict = defaultdict(lambda:defaultdict(int)) #dict mapping right hand side of rule to left hand sides  or rule to integer frequency
    rhsFreq = defaultdict(int) #dict mapping rhs sides to total frequencies
    bDict = defaultdict(lambda:defaultdict(int)) #dict mapping Bs to A,B,C touples to frequencies for binary rules of the type A -> BC
    uniDict = defaultdict(lambda:defaultdict(int)) #maps uniary rule to its frequency
    uniFreq = defaultdict(int) #frequency fo uniary rules 

    for line in data.readlines():
        if len(line.split()) == 5: #b if the rule is a binary  rule
            count = int(line.split()[0])
            a = str(line.split()[1])
            b = str(line.split()[3])
            c = str(line.split()[4])
            rhside = str(line.split()[3]) + " " + str(line.split()[4])
            tripple = a, b, c
            bDict[b][tripple] = count
            count = int(line.split()[0])
            rhsFreq[rhside] += count
            rhsDict[rhside][a] = count
        if len(line.split()) == 4: #uniary rule
            count = int(line.split()[0])
            bside = line.split()[3] #right hand side of rule
            aside = line.split()[1]# left hand side of rule
            uniDict[bside][aside] = count
            uniFreq[bside] += count
    for bsides in uniDict.keys():
        for asides in uniDict[bsides].keys():
            if uniDict[bsides][asides] == 1: #delete all the uniary rules that only occour once
                del uniDict[bsides][asides]
                uniFreq[bsides] = uniFreq[bsides] - 1
    print "done buiding dicts"

    return (bDict, rhsDict, rhsFreq,uniDict, uniFreq)

class Constituent:
    def __init__(self, label, pointer1, pointer2, mew):
        self.label = label
        self.pointer1 = pointer1
        self.pointer2 = pointer2
        self.mew = mew
class Cell:
    def __init__(self):
        self.constitList = []
        self.cellDict = { } #blank list of topics to get over written each time. Used to remove old topics

def addUniary (C,i,k, rhsDict, rhsFreq, uniDict, uniFreq): #creates consituents for all possible uniary rules
    #print "starting uniary additions for cell", i, k
    recently_added  = C[i][k].constitList
    prev_added = []

    while len(recently_added) > 0: #keep adding consituents from uniary rules until an iteration generates no new constituents
        temp = []
        for constit in recently_added: #only consitutents from the most recent iteration have the potential to generate new constituents
            count = float(uniFreq[constit.label])
            bside = constit.label
            for aside in uniDict[bside].keys(): #add all possible As
                tempProb =  float(uniDict[bside][aside]) / count) 
                thisProb = math.log(tempProb) + constit.mew # Pr * mewC1  * mewC2

                for ogConstit in C[i][k].constitList: #loop thru the original labels for potential A sides
                    if ogConstit.label == aside and ogConstit.mew < thisProb:
                        C[i][k].constitList.remove(ogConstit) #impossible for this consit to be the best constit w/ this label
                    if ogConstit.label == aside and ogConstit.mew >= thisProb: #already a better consist w/ this label in the cell
                        break
                else: #a better consist w/ this label doesn't exist, so we need to create it
                    newConstit = Constituent(aside,constit,None,thisProb)
                    C[i][k].constitList.append(newConstit)
                    C[i][k].cellDict[aside] = newConstit
                    temp.append(newConstit)
                    #print "uniary", newConstit.label, "and prob", newConstit.mew, "to C", i, k
        prev_added += recently_added
        recently_added = temp #overwrite temp w/ consits added at this iteration 

def fill (C,i,k,string, bDict, rhsDict, rhsFreq, uniDict, uniFreq, L):
    #print "starting to work on filling span", i, k
    if k == i + 1: # w(i:k) is a one word terminal
        word = str(string.split()[i])
        newConstit = Constituent(word, None, None, 0.000) #terminal constit, log(1) = 0
        C[i][k].constitList.append(newConstit)
        C[i][k].cellDict[word] = newConstit #add the terminal word 
        addUniary(C,i,k, rhsDict, rhsFreq, uniDict, uniFreq) #add new consists generated from uniary rules off of terminal word
    else: #not a terminal
        for j in range (i +1, k):
            for leftChild in C[i][j].constitList:
                #print "looking at left child wtih label", leftChild.label
                for tripple in bDict[leftChild.label]: #for rules of a-->bc, look at all cs for this b
                    aside = tripple[0]
                    bside = tripple[1]
                    cside = tripple[2]
                    bc = bside + " " + cside
                    if C[j][k].cellDict.has_key(cside):
                        rightChild = C[j][k].cellDict[cside]
                        count = float(rhsFreq[bc]) # number of times any a generates this bc 
                        newProb =  math.log( float(rhsDict[bc][aside]) / count) + leftChild.mew + rightChild.mew
                        if len(C[i][k].constitList) == 0: #if the cell we are filling is empty
                            newConstit = Constituent(aside, leftChild, rightChild, newProb)
                            C[i][k].constitList.append(newConstit)
                            C[i][k].cellDict[aside] = newConstit
                            #print "added initial binary consit w/ label", aside, "p vec = ", newProb, "at", i, k
                        else:  # this cell is not empty 
                            for priorConstit in C[i][k].constitList: # go over the constits already in the cell
                                if priorConstit.label == aside and priorConstit.mew < newProb: # can't be the best constit with this label
                                    C[i][k].constitList.remove(priorConstit)
                                if priorConstit.label == aside and priorConstit.mew > newProb: #already a better consist w/ that label
                                    break
                            else:
                                newConstit = Constituent(aside, leftChild, rightChild, newProb)
                                C[i][k].constitList.append(newConstit)
                                C[i][k].cellDict[aside] = newConstit
                                #print "ADDED CONSIT AT", i, k, "WITH P VEC", newProb, "AND LABEL", newConstit.label

        addUniary(C,i,k, rhsDict, rhsFreq, uniDict, uniFreq)


def parse(string,C, bDict, rhsDict,rhsFreq, uniDict, uniFreq, L):
    toParse = []
    for word in string.split():
        toParse.append(word)
    L = len(toParse)
    for l in range(1, L+1):
        for s in range (0, ((L-l)+1)):
            fill(C,s,s+l,string, bDict, rhsDict, rhsFreq, uniDict, uniFreq, L)



def trace(constituent, visited): # prefrom a BFS on cell pyriamd to read back out the sentence tree digram
    visited += "("
    visited += constituent.label
    if constituent.pointer1 != None:
        visited = trace(constituent.pointer1, visited)
    if constituent.pointer2 != None:
        visited = trace(constituent.pointer2, visited)
    visited += ")"
    return visited





if __name__ == "__main__":
    file1 = sys.argv[1] # rule blan k
    toParse = sys.argv[2] #wsj-24.txt
    output = sys.argv[3]
    bDict,rhsDict, rhsFreq, uniDict,uniFreq = buildDicts(file1)
    toParse = open(toParse, "r")
    output = open(output, 'w')
    sentencenumber = 0

    for sentence in toParse.read().split(' .'):
        sentencenumber += 1
        print sentencenumber

        count = 0
        for word in sentence.split():
            count += 1
        if count >= 25:
            output.write("*IGNORE*" + "\n")
        if count < 25:
            sentence = sentence + " " + "."
            #print sentence, "what i think the sentence im parsing is"
            C = [[] for i in range(count + 2)]
            for row in C:
                for i in range(count+2): #probably a more elegant way to do this #create a list of empty lists the size of the length of the sentence we are trying to parse 
                    newCell = Cell()
                    row.append(newCell)
            parse(sentence, C, bDict, rhsDict, rhsFreq,uniDict, uniFreq, count)

            bestProb = -999999999
            bestConstit = Constituent("*UNK*",None,None, -9999999)
         
            if (C[0][count +1].constitList) == None:
                output.write("*IGNORE*")
                print "AN ERROR HAS OCCURED"
            else:

                for constit in C[0][count + 1].constitList:  #+1 becuase period?
                    if constit.mew > bestProb and constit.label == "TOP":
                        bestConstit = constit
                        bestProb = constit.mew

                        visited = trace(bestConstit, "")
                output.write("\n" + visited)
                print "\n" + visited



    