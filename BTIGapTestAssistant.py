# Author: Kylie Hoyt
# Date Created: 06/03/2022
# How to use this script: This script works in conjunction with 3 .csv files: a varaible database, a test plan, and a sample table. The variable database contains the part numbers of all the equipment to be randomized including the divoted implants, posterior cutting blocks, and the drill through trials. The variable database also contains the workflow options to be randomized and the remaining uses of each of the burs being used. The test plan contains the approximate number of Navio and CORI cases performed by each user and the samples that each user still has to perform. The plan for this test is 4 saw samples, 4 5mm cylindrical bur samples, and 4 5mm bullet bur samples per user. The sample table contains the conditions of each sample including the user information, the randomized variables, and the serial numbers of the CORI equipment being used. These tables are uploaded when the script begins and updated/saved after every confirmed sample and at the end. The real time randomization takes into account the innaccessibility of each divoted implant while it is waiting for divot collection or pull-off testing. This feature can be disabled and a randomization matrix can be generated instead.
# This code was created for Smith & Nephew to help run a feasibility test called bone-to-implant gap testing. The test administrator should use this code to maintain the sample information and keep track of available equipment. 
import numpy as np
from numpy.random import randint
from numpy import random
import pandas as pd

# Serial Numbers need to be manually updated if they change
ConsoleSN = 2231
CameraSN = 4124
HandpieceSN = 2543
LongAttachmentSN = 6536
BoneHardwareIdentifier = 4635
BoneHardwareAN = 32432
FootpedalSN = 352352
TabletSN = 3525
DisplaySN = 45745
CartSN = 63463
BurMethods = {'6mm Bullet Bur': 0, '6mm Cylindrical Bur': 1} # Maps the names of the prep methods to the indices of the methods for tracking bur uses in the Variable Database
TibiaWorkFlowMessages = {'Hybrid': 'Hybrid: \n1) Use the 3mm spacer to place the tibia \n2) Bur the distal and posterior chamfer surfaces \n3) Either bur the tibia or use the 12mm spacer to simulate burring the tibia \n4) Bur the posterior surface', 'Femur First': 'Femur First: \n1) Use the 3mm spacer to place the tibia \n2) Bur all surfaces of the femur', 'Tibia First Bur': 'Tibia First: \n1) Either use the 3mm spacer to place the tibia then bur the tibia or use the 12mm spacer to simulate burring the tibia \n2) Bur all surfaces of the femur', 'Tibia First Saw': 'Engage Method: \n1) Either use the 3mm spacer to place the tibia then saw the tibia or use the 12mm spacer to simulate sawing the tibia \n2) Saw all surfaces of the femur'}


class User:
    # def __init__ (self, Username, Samples, CaseRecord):
    # Inputs: 
    #   Username - Name of the user as defined in the Test Plan
    #   Samples - The names of the remaining samples the user needs to perform
    #   CaseRecord - [Approximate number of Navio case performed, Approximate number of CORI cases performed]
    def __init__(self, Username, Samples, CaseRecord):
        self.Name = Username
        self.RemainingSamples = Samples
        self.NumNavioCases = CaseRecord[0]
        self.NumCoriCases = CaseRecord[1]
    # def GenerateSample(self):
    # Purpose: Randomly generate a sample prep plan and update the sample table, test plan, and equipment list accordingly
    def GenerateSample(self):
        # Randomly choose a sample prep method from the user's remaining samples
        if len(self.RemainingSamples) == 0:
            print("This user has no remaining samples.")
            return
        RandomSampleIndex = randint(len(self.RemainingSamples))
        PrepMethod = self.RemainingSamples[RandomSampleIndex]
        # Randomly choose a divoted implant from the available implants
        if len(Variables.DivotedImplants) == 0:
            print("No divoted implants available.")
            Variables.CheckDivotedImplantInventory()
            if len(Variables.DivotedImplants) == 0:
                print("No divoted implants available. No sample can be generated")
                return
        DivotedImplantIndex = randint(len(Variables.DivotedImplants))
        DivotedImplantPN = Variables.DivotedImplants[DivotedImplantIndex]
        # If the sample is to be prepped with a saw, a posterior cutting block must be randomly chosen
        if PrepMethod == 'Saw':
            TibiaWorkflow = Variables.SawWorkflows[randint(len(Variables.SawWorkflows))]   
            # Randomly choose a posterior cutting block
            PosteriorCuttingBlockPN = Variables.PosteriorCuttingBlocks[randint(len(Variables.PosteriorCuttingBlocks))]
            NewSample = Sample('Saw', DivotedImplantPN, TibiaWorkflow, "N/A", PosteriorCuttingBlockPN)
            # Confirm sample plan
            NewSample.PrintSampleSummary()
            print("Is necessary equipment available? (y/n)")
            Answer = input()
            # If user wants to proceed with sample
            if Answer.lower() == 'y':
                self.RemainingSamples.pop(RandomSampleIndex)    # Remove sample from user's remaining sample list
                SampleUser = User(self.Name, self.RemainingSamples, [self.NumNavioCases, self.NumCoriCases])
                CurrentSamples.AddSample(SampleUser, NewSample) # Add sample to the sample table
                Variables.MakeDivotedImplantInaccessible(DivotedImplantIndex)
                BTIGapTest.UpdatePlan(SampleUser)               # Update the user in the test plan
                Variables.SaveVariableDatabaseAsCSV()
                BTIGapTest.SavePlanAsCSV()
                CurrentSamples.SaveSampleTableAsCSV()
                print("Data Saved\n\n")
        # If the sample is to be prepped with a bur, a tibia placement and drill through trial must be randomly chosen, and bur use counts must be updated
        if PrepMethod == '6mm Bullet Bur' or PrepMethod == '6mm Cylindrical Bur':
            # Randomly choose a workflow to guide the placement of the tibia
            TibiaWorkflow = Variables.BurWorkflows[randint(len(Variables.BurWorkflows))]
            # Randomly choose a drill through trial
            DrillThroughTrialPN = Variables.DrillThroughTrials[randint(len(Variables.DrillThroughTrials))]
            NewSample = Sample(PrepMethod, DivotedImplantPN, TibiaWorkflow, DrillThroughTrialPN, "N/A")
            # Confirm sample plan
            NewSample.PrintSampleSummary()
            Variables.CheckBurUses(PrepMethod)  # Alert if bur has been used 4 times
            print("Is necessary equipment available? (y/n)")
            Answer = input()
            if Answer.lower() == 'y': 
                self.RemainingSamples.pop(RandomSampleIndex)    # Remove sample from user's remaining sample list
                SampleUser = User(self.Name, self.RemainingSamples, [self.NumNavioCases, self.NumCoriCases])
                CurrentSamples.AddSample(SampleUser, NewSample) # Add sample to sample table
                Variables.MakeDivotedImplantInaccessible(DivotedImplantIndex)
                Variables.UseBur(PrepMethod)                    # Keep track of bur uses
                BTIGapTest.UpdatePlan(SampleUser)               # Update the user in the test plan
                Variables.SaveVariableDatabaseAsCSV()
                BTIGapTest.SavePlanAsCSV()
                CurrentSamples.SaveSampleTableAsCSV()
                print("Data Saved\n\n")


class Sample:
    # def __init__(self, Method, ImplantPN, Workflow, TrialPN, CutBlockPN):
    # Inputs:
    #   Method - Name of the prep method being used
    #   ImplantPN - Part number of the divoted implant being impacted
    #   Workflow - Key of the workflow being used to guide the placement of the tibia
    #   TrialPN - Part numnber of the drill through trial being used; N/A for saw-prepped samples
    #   CutBlockPN - Part number of the posterior cutting block being used; N/A for bur-prepped samples
    def __init__(self, Method, ImplantPN, Workflow, TrialPN, CutBlockPN):
        self.PrepMethod = Method
        self.DivotedImplantPN = ImplantPN
        self.TibiaWorkflow = Workflow
        self.DrillThroughTrialPN = TrialPN
        self.PosteriorCuttingBlockPN = CutBlockPN
    def PrintSampleSummary(self):
        print("\n\nPrep Method: ", self.PrepMethod)
        print("Workflow: ", TibiaWorkFlowMessages[self.TibiaWorkflow])
        if self.PrepMethod == 'Saw':
            print("Posterior Cutting Block PN: ", self.PosteriorCuttingBlockPN)
        elif self.PrepMethod == '6mm Cylindrical Bur' or self.PrepMethod == '6mm Bullet Bur':
            print("Drill Through Trial PN: ", self.DrillThroughTrialPN)
        print("Divoted Implant PN: ", self.DivotedImplantPN)

        


class VariableDatabase:
    # def __init__(self, path):
    # Inputs:
    #   path - full file path of the location of the variable database csv file; Column names must match the names below
    def __init__(self, path):
        self.filepath = path
        Data = pd.read_csv(self.filepath)        
        self.DivotedImplants = list(Data["DivotedImplantPNs"].dropna().astype(int).astype(str))
        self.PosteriorCuttingBlocks = list(Data["PosteriorCuttingBlockPNs"].dropna().astype(int).astype(str))
        self.DrillThroughTrials = list(Data["DrillThroughTrialPNs"].dropna().astype(int).astype(str))
        self.BurWorkflows = Data["BurWorkflows"].dropna()
        self.SawWorkflows = Data["SawWorkflows"].dropna()
        self.BurUseCounts = pd.to_numeric(Data["BurRemainingUses"].dropna())
    def MakeDivotedImplantInaccessible(self, ImplantIndex):
        self.DivotedImplants.pop(ImplantIndex)
        print("Remaining Divoted Implants: ", self.DivotedImplants)
    def MakeDivotedImplantAccessible(self, ImplantPN):
        self.DivotedImplants = self.DivotedImplants + [ImplantPN]
    def CheckDivotedImplantInventory(self):
        while True:
            print("Need to add available divoted implant? (y/n):")
            AddDivotedImplant = input()
            if AddDivotedImplant == 'y':
                print("Enter divoted implant PN: ")
                ImplantPN = input()
                self.MakeDivotedImplantAccessible(ImplantPN)
            else:
                return
    def UseBur(self, PrepMethod):
        # BurMethods[PrepMethod] maps the name of the method to the index of the method in the variable database
        self.BurUseCounts[BurMethods[PrepMethod]] = self.BurUseCounts[BurMethods[PrepMethod]] - 1
    def CheckBurUses(self, PrepMethod):
        if self.BurUseCounts[BurMethods[PrepMethod]] == 0:
            print("**********Bur has been used 4 times. Replace bur before proceeding**********")
            self.BurUseCounts[BurMethods[PrepMethod]] = 4
    def Summarize(self):
        print("Divoted Implant PNs: \n", self.DivotedImplants)
        print("Drill Through Trial PNs: \n", self.DrillThroughTrials)
        print("Posterior Cutting Block PNs: \n", self.PosteriorCuttingBlocks)
    def SaveVariableDatabaseAsCSV(self):
        VariableDict = {"DivotedImplantPNs": self.DivotedImplants, "PosteriorCuttingBlockPNs": self.PosteriorCuttingBlocks, "DrillThroughTrialPNs": self.DrillThroughTrials, "BurWorkflows": self.BurWorkflows, "SawWorkflows": self.SawWorkflows, "BurRemainingUses": self.BurUseCounts}
        df = pd.DataFrame(dict([(key,pd.Series(val)) for key,val in VariableDict.items()]))
        df.to_csv(self.filepath, index=False)


class TestPlan:
    # def __init__(self, path):
    # Inputs:
    #   path - full file path of the location of the test plan csv file; Column names are the names of the users
    def __init__(self, path):
        self.filepath = path
        Data = pd.read_csv(self.filepath, dtype=str)
        self.Users = {}
        # In each column, the first entry is the number of Navio cases performed, the second entry is the number of CORI cases performed, and the remaining entries are the samples the user still has to perform
        for column, value in Data.iteritems():
            self.Users[column] = list(Data[column].dropna())
            NewUser = User(column, self.Users[column][2:], self.Users[column][0:2]) # User(Name, RemainingSamples, CaseRecord)
    def UpdatePlan(self, UpdatedUser):
        UserInfo = [UpdatedUser.NumNavioCases, UpdatedUser.NumCoriCases] + list(UpdatedUser.RemainingSamples)
        self.Users[UpdatedUser.Name] = UserInfo
        print("Remaining Samples for ", UpdatedUser.Name, ": ", UpdatedUser.RemainingSamples)
    def SummarizePlan(self):
        print("Remaining Samples in Test Plan:")
        for key, val in self.Users.items(): # User Name: list of RemainingSamples
            print(key, ": \n", val[2:])
        return
    def SavePlanAsCSV(self):
        df = pd.DataFrame({ key: pd.Series(val) for key, val in self.Users.items()})
        df.to_csv(self.filepath, index=False)


class SampleDatabase:
    # def __init__(self, path):
    # Inputs:
    #   path - full file path of the location of the sample table csv file; Column names are the number of the sample; The first column is the names of the rows
    def __init__(self, path):
        self.filepath = path
        self.Samples = pd.read_csv(self.filepath, na_filter=False)
    def AddSample(self, SampleUser, SampleInfo):
        # Append a new column to the self.Samples dataframe
        SampleNumber = len(self.Samples.columns)
        self.Samples[SampleNumber] = [SampleUser.Name, SampleUser.NumNavioCases, SampleUser.NumCoriCases, SampleInfo.PrepMethod, SampleInfo.DivotedImplantPN, SampleInfo.TibiaWorkflow, SampleInfo.DrillThroughTrialPN, SampleInfo.PosteriorCuttingBlockPN, ConsoleSN, CameraSN, HandpieceSN, LongAttachmentSN, BoneHardwareIdentifier, BoneHardwareAN, FootpedalSN, TabletSN, DisplaySN, CartSN]
        print("\n\nNewest Sample:")
        print(self.Samples[['Sample Number', SampleNumber]]) # Show names of the rows and the latest sample
    def GenerateRandomTable(self):
        SampleNumber = 0
        TempSamples = pd.DataFrame()
        for us, val in BTIGapTest.Users.items():
            for v in val[2:]: # for each prep method in remaining samples
                SampleNumber += 1
                DivotedImplantIndex = randint(len(Variables.DivotedImplants))
                DivotedImplantPN = Variables.DivotedImplants[DivotedImplantIndex]
                # If the sample is to be prepped with a saw, a posterior cutting block must be randomly chosen
                if v == 'Saw':
                    TibiaWorkflow = Variables.SawWorkflows[randint(len(Variables.SawWorkflows))]   
                    # Randomly choose a posterior cutting block
                    PosteriorCuttingBlockPN = Variables.PosteriorCuttingBlocks[randint(len(Variables.PosteriorCuttingBlocks))]
                    DrillThroughTrialPN = "N/A"
                elif v == '6mm Bullet Bur' or v == '6mm Cylindrical Bur':
                    # Randomly choose a workflow to guide the placement of the tibia
                    TibiaWorkflow = Variables.BurWorkflows[randint(len(Variables.BurWorkflows))]
                    # Randomly choose a drill through trial
                    DrillThroughTrialPN = Variables.DrillThroughTrials[randint(len(Variables.DrillThroughTrials))]
                    PosteriorCuttingBlockPN = "N/A"
                TempSamples[SampleNumber] = [us, val[0], val[1], v, DivotedImplantPN, TibiaWorkflow, DrillThroughTrialPN, PosteriorCuttingBlockPN, ConsoleSN, CameraSN, HandpieceSN, LongAttachmentSN, BoneHardwareIdentifier, BoneHardwareAN, FootpedalSN, TabletSN, DisplaySN, CartSN]
        SortedDf = TempSamples.sort_values(by=[3, 0], axis=1) # sort columns by prep method then user
        LabeledDf = pd.concat([self.Samples, SortedDf], axis=1) # add columns to existing sample table
        LabeledDf.to_csv(self.filepath, index=False, na_rep='N/A')
        print(LabeledDf)
    def SaveSampleTableAsCSV(self):
        self.Samples.to_csv(self.filepath, index=False, na_rep='N/A')


# Upon Running

#Start of Day
print("Enter filepath for Variable Database:")
Variablefile = input()
Variables = VariableDatabase(Variablefile)
Variables.Summarize()
print("Enter filepath for Test Plan:")
TestPlanfile = input()
BTIGapTest = TestPlan(TestPlanfile)
BTIGapTest.SummarizePlan()
print("Enter filepath for Sample Table:")
SampleTablefile = input()
CurrentSamples = SampleDatabase(SampleTablefile)
print(CurrentSamples.Samples, "\n\n")

##########Matrix Mode##########
CurrentSamples.GenerateRandomTable()

##########Interactive Mode##########
# Variables.CheckDivotedImplantInventory()
# print("Run a new sample? (y/n): ")
# KeepRunning = input().lower()
# while KeepRunning == 'y':
    # print("Enter User Name:")
    # Name = input()
    # SampleUser = User(Name, BTIGapTest.Users[Name][2:], BTIGapTest.Users[Name][0:2])
    # SampleUser.GenerateSample()
    # Variables.CheckDivotedImplantInventory()
    # print("Run a new sample? (y/n): ")
    # KeepRunning = input().lower()
# print("Current Sample Table:")
# print(CurrentSamples.Samples)
# CurrentSamples.SaveSampleTableAsCSV()	
# print("\n\nSample Collection Stopped")
    
#End of Day
BTIGapTest.SummarizePlan()
Variables.SaveVariableDatabaseAsCSV()
BTIGapTest.SavePlanAsCSV()
print("Data Saved")
