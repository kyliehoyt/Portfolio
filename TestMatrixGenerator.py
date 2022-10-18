# Author: Kylie Hoyt
# Date Created: 06/27/2022
# How to use this script:
# This script is used to randomly generate a sample table given a sample
# space and values for variables required to track. The script works in
# conjunction with 2 .txt files: a sample space file and a variable
# file. The sample space file contains the Users, independent variables,
# desire for random trial order, sample groups, sample group replicates,
# and the variables required by each sample group. The variables file
# contains the random variables, constants, fields, and their column
# names for the sample table. The script reads in both files and
# populates a sample table with random variables, constants, and fields
# specific to each sample group for each user and sample group
# replicate.
import numpy as np
from numpy.random import randint
import random
import pandas as pd
import os
import csv
import ast
import sys


def pathStrip(input,levels=0): # cd .. from file path input for as many levels as specified
        output = input
        for i in range (levels + 1):
            output, tmp2 = os.path.split(output)
        return output


def randomOrder(sortedSampleTable):
        # create a random order
        newIndex = list(range(len(sortedSampleTable)))
        random.shuffle(newIndex)
        shuffledIndexDataframe = pd.DataFrame({"newIndex": newIndex}, index=None)
        labeledDataframe = pd.concat([shuffledIndexDataframe, sortedSampleTable], axis=1)
        # sort samples by new random order
        randomDataframe = labeledDataframe.sort_values(by=["User", "newIndex"], axis=0, ignore_index=True)
        # relabel samples with ordered sample number
        randomDataframe.drop(["newIndex", "Sample"], axis=1, inplace=True)
        sampleIndex = list(range(1, len(sortedSampleTable)+1))
        sortedIndexDataframe = pd.DataFrame({"Sample": sampleIndex}, index=None)
        return pd.concat([sortedIndexDataframe, randomDataframe], axis=1)

        

if len(sys.argv)==4:
        sampleSpaceFilePath = sys.argv[1]
        variablesFilePath = sys.argv[2]
        sampleTableFileName = sys.argv[3]
else:
        print("Enter filepath for Sample Space:")
        sampleSpaceFilePath = input()
        print("Enter filepath for Variable Database:")
        variablesFilePath = input()
        print("Enter Desired Sample Table File Name:")
        sampleTableFileName = input()


sampleSpaceFile = open(sampleSpaceFilePath, 'r')
sampleSpaceContent = sampleSpaceFile.readlines(-1)
sampleSpace = {} # dictionary containing the users, independent variables, and variables required for each sample group
for line in sampleSpaceContent:
        if "#" in line:
                continue
        head,sep,tail = line.partition('=')
        sampleSpace[head.strip()] = ast.literal_eval(tail.strip('\n').strip())
sampleSpaceFile.close()


variablesFile = open(variablesFilePath, 'r')
variablesContent = variablesFile.readlines(-1)
variables = {} # dictionary containing the variables, their possible values, sample table column names, and fields
for line in variablesContent:
        if "#" in line:
                continue
        head,sep,tail = line.partition('=')
        variables[head.strip()] = ast.literal_eval(tail.strip('\n').strip())
variablesFile.close()


sampleTableFile = pathStrip(variablesFilePath)
sampleTableFilePath = sampleTableFile + "\\" + sampleTableFileName + ".csv"

# Making the table
sampleGroups = list(sampleSpace.keys())[3:] # first three keys are Users, Independent Variables, Random Order
sampleColumnNames = ['Sample', 'User'] + sampleSpace["Independent Variables"]
# column names are the first elem of each variable list; different variables may be presented in the same column
varColumnNames = []
[varColumnNames.append(elem[0]) for elem in list(variables.values()) if elem[0] not in varColumnNames]
colList = sampleColumnNames + varColumnNames
# empty sample table
tempDict = {}
data = pd.DataFrame(tempDict, columns=colList)


sample = 0
for user in sampleSpace['Users']: # for each user
        for method in sampleGroups: # for each combination of independent variables
                for n in range(sampleSpace[method][0]): # for the number of replicates for each sample group
                        sample = sample+1
                        # start df row with sample number, user, and sample group
                        sampleInfo = [sample, user] + method.split("-")
                        sampleData = {k:v for (k,v) in zip(sampleColumnNames, sampleInfo)}
                        for var,vals in variables.items(): # for each random variable, constant, or field
                                if var in sampleSpace[method]: # if the sample group requires it
                                        if len(vals) > 1: # random variable or constant needed for sample group
                                                valIndex = randint(1,(len(vals))) # randomly choose a value
                                                val = vals[valIndex]
                                        else: # field needed for sample group
                                                val = ' '
                                elif vals[0] in sampleData.keys(): # if it's not required but the column already exists
                                        continue # don't overwrite the shared column
                                else: # variable not needed for sample group
                                        val = 'N/A'
                                sampleData[vals[0]] = val # add column to the row
                        # append row to bottom of sample table
                        sampleData = pd.DataFrame(sampleData, index=[sample-1], columns=colList) 
                        data = pd.concat([data, sampleData], axis=0)

if sampleSpace["Random Order"] == True:
        data = randomOrder(data)

data.to_csv(sampleTableFilePath, index=False, na_rep='N/A')
