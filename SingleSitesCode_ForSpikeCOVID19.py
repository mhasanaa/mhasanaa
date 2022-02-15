# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 17:00:28 2021

@author: Hasan
"""

import xlrd
import pandas as pd
import time
import numpy as np
import statistics as st
import math
from datetime import datetime
from openpyxl import load_workbook

start_time = time.process_time()
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Start =", current_time)

def count_dups(nums):
    element = []
    freque = []
    nums = sorted(nums)
    if not nums:
        return element
    running_count = 1
    for i in range(len(nums)-1):
        if nums[i] == nums[i+1]:
            running_count += 1
        else:
            freque.append(running_count)
            element.append(nums[i])
            running_count = 1
    freque.append(running_count)
    element.append(nums[i+1])
    return element,freque

def concat_dups(criteria,ele):
    element = []
    freque = []
    if not criteria:
        return element
    
    running_count = ele[0]
    for i in range(len(criteria)-1):
        if criteria[i] == criteria[i+1]:
            running_count += ","+ele[i+1] #the original form is without str
        else:
            freque.append(running_count)
            element.append(criteria[i])
            running_count = ele[i+1]

    freque.append(running_count)
    element.append(criteria[i+1])
    return element,freque

def sumsimilar(criteria,ele):
    element = []
    freque = []
    if not criteria:
        return element
    
    running_count = ele[0]
    for i in range(len(criteria)-1):
        if criteria[i] == criteria[i+1]:
            running_count += ele[i+1]
        else:
            freque.append(running_count)
            element.append(criteria[i])
            running_count = ele[i+1]

    freque.append(running_count)
    element.append(criteria[i+1])
    return element,freque

def updateExcel(df_append, file_name):
    to_update = {"Sheet1": df_append} # define what sheets to update
    excel_reader = pd.ExcelFile(file_name) # load existing data
    excel_writer = pd.ExcelWriter(file_name)# write and update

    for sheet in excel_reader.sheet_names:
        sheet_df = excel_reader.parse(sheet)
        append_df = to_update.get(sheet)

        if append_df is not None:
            sheet_df = pd.concat([sheet_df, append_df], axis=1)
        
        sheet_df.to_excel(excel_writer, sheet, index=False)

    excel_writer.save()

def strip_n(f1):
    f2=[]
    for i in range(len(f1)):
        f1[i]=f1[i].rstrip('\r\n')
        f2.append(f1[i])
    return f2

def raw_WL(ref0,ProtLen,seq):

    OriAA, Pos, Residues, MutAA, Count, WL_Count, NC= [],[],[],[],[],[],[]

    for i in range(ProtLen):#Loop every position
        OriAA_Pass, Pos_Pass, Residues_Pass, MutAA_Pass, Count_Pass, WL_Pass, NC_Pass = [],[],[],[],[],[],[]  #to store AA for each seq
        list_AA = []
        # ref0 = f0[0] # store original AA
    
        for j in range(len(seq)):#Loop every seq
            seq0 = seq[j] #assign a seq
            aa0 = seq0[i] #record the AA of the seq on the corresponding position
            list_AA.append(aa0)
    
        AA, count = count_dups(list_AA)
        for k in range(len(AA)):
            if AA[k] != ref0[i]: #Comment this line if you want to generate plot for full version
                MutAA_Pass.append(AA[k])# list of mutated position
                Count_Pass.append(count[k])#list of number of seq from mutated pos
        
        for k in range(len(MutAA_Pass)):#Record the original residues
            add = ref0[i]+str(i+1)
            OriAA_Pass.append(ref0[i])
            Residues_Pass.append(add)#list of original residues
            Pos_Pass.append(int(i+1))# list of position
            WL_Pass.append(MutAA_Pass[k]+"("+str(Count_Pass[k])+")")
            NC_Pass.append(add + MutAA_Pass[k])
    
        if len(MutAA_Pass) != 0:
            OriAA.extend(OriAA_Pass)#Original Amino Acid
            Pos.extend(Pos_Pass)#Position of Mutated AA
            Residues.extend(Residues_Pass)#Position + Original AA
            MutAA.extend(MutAA_Pass)#Mutated AA for Each Position 
            Count.extend(Count_Pass)#Mutation Count for each MutAA
            WL_Count.extend(WL_Pass)#The Weblogo/single Site format
            NC.extend(NC_Pass)#Mutation format
    
    return OriAA,Pos,Residues,MutAA,Count,WL_Count,NC #OriAA = Ori

def WL(OriAA,Pos,Residues,MutAA,Count,WL_Count,NC,ref0,dfTypeAA):
    d = {'Ori':OriAA,'Pos': Pos,'Residues': Residues,'MutAA': MutAA,'Mutation': NC, \
         'Count': Count,'WL-Count':WL_Count,};  
    df = pd.DataFrame(d)

    df = pd.merge(df, dfTypeAA, on ='Mutation',how ='left')#Merging the information about AA Type
    df["NC"] = df["NC"].astype(str)
    df["space"] = "" #df = The raw data
    df = df.fillna(0)
    
    FullResidues = []
    for m in range(len(ref0)):
        w = ref0[m]+str(m+1)
        FullResidues.append(str(w))
    
    # print(FullResidues)

    Residue= df["Residues"].tolist()
    MutAA = df["MutAA"].tolist()
    Residue, MutAA = concat_dups(Residue,MutAA) #WL_AA

    Residue2= df["Residues"].tolist()
    WL_Count = df["WL-Count"].tolist()
    Residue2, WL_Count = concat_dups(Residue2,WL_Count) #WL_count

    Residues = df["Residues"].tolist()
    Residues, ResCount = count_dups(Residues) #N.AA

    Residues2 = df["Residues"].tolist()
    NCount = df["Count"].tolist()
    Residues2, NCount = sumsimilar(Residues2,NCount)#N.Count
    

    d1 = {'Sites': FullResidues}; df1 = pd.DataFrame(d1)
    d2 = {'Sites': Residues,'N.AA': ResCount}; df2 = pd.DataFrame(d2)
    df3 = pd.merge(df1, df2, on ='Sites',how ='left')
    
    
    d30 = {'Sites': Residues2,'N.Count': NCount}; df30 = pd.DataFrame(d30)
    df3 = pd.merge(df3, df30, on ='Sites',how ='left')
    
    d31 = {'Sites': Residue,'WL-AA': MutAA}; df31 = pd.DataFrame(d31)
    df3 = pd.merge(df3, df31, on ='Sites',how ='left')
    
    d32 = {'Sites': Residue2,'WL-Count': WL_Count}; df32 = pd.DataFrame(d32)
    df3 = pd.merge(df3, df32, on ='Sites',how ='left')
    
    df3 = df3.fillna(0)
     
    
    dfNAA = df3[['Sites', 'N.AA']].copy()
    dfWL_AA = df3[['Sites', 'WL-AA']].copy()
    dfNCount = df3[['Sites', 'N.Count']].copy()
    dfWL_Count = df3[['Sites', 'WL-Count']].copy()
    # print(dfNAA)
    # print(dfWL_AA)
    # print(dfNCount)
    
    
    df3["Pos"] =list(range(1,ProtLen+1))
    df3["space"] = ""
    
    return df,df3,dfNAA,dfNCount,dfWL_AA,dfWL_Count #,dfNCount # df = Raw data of WL ; df3 = Final raw data

def ProbDist(OriAA,Residues,MutAA,Count,df1,RefGene):# Average for each sites
    d = {'Ori':OriAA,'Residues': Residues, \
         'MutAA': MutAA,'Count': Count};  df = pd.DataFrame(d)
    df = pd.merge(df, df1, on ='Residues',how ='left') #df1 contain residues and Ncount column
    
    Count = df["Count"].tolist(); NCount = df["N.Count"].tolist()
    Ori = df["Ori"].tolist(); Mut = df["MutAA"].tolist()
    Residues = df["Residues"].tolist()
    
    Prob,OriMut = [],[]
    FinProb,Residues2 = [],[]
    for i in range(len(Count)):
        Prob.append(round(Count[i]/NCount[i],5))
        OriMut.append(Ori[i]+Mut[i])# OriMut ==> to identify different type permutation
    
    df["OriMut"] = OriMut;df["Prob"] = Prob
    df2 = df.filter(['OriMut','Prob'], axis=1)
    df2 = df2.groupby(['OriMut']).sum()#Numerator
    df2 = df2.reset_index(level=['OriMut']) #To revive the previous summed up OriMut
    
    OriMut = df2["OriMut"].tolist()
    Residues2.extend([i[:1] for i in OriMut])
    df2["Residues"] = Residues2
    
    
    df3 = df2.filter(['Residues','Prob'], axis=1)
    df3 = df3.groupby(['Residues']).sum()#Denumerator
    df3 = df3.reset_index(level=['Residues']) #To revive the previous summed up Residues
    # 
    
    df2 = pd.merge(df2, df3, on ='Residues',how ='left')#final probability
    Num = df2["Prob_x"].tolist();Den = df2["Prob_y"].tolist();Permut = df2["OriMut"].tolist()

    for i in range(len(Num)):
        FinProb.append(round(Num[i]/Den[i],5))
    
    d4 = {'Permut':Permut,'Prob': FinProb};  df4 = pd.DataFrame(d4)
    d5 = {'Permut':RefGene};  df5 = pd.DataFrame(d5)
    df5 = pd.merge(df5, df4, on ='Permut',how ='left')#final data frame
    df5 =df5.fillna(0)
    print(df5)
    
    Prob = df5["Prob"].tolist()
    x = int(np.sqrt(len(RefGene))); y = int(np.sqrt(len(RefGene)))
    Prob = np.array(Prob).reshape((x,y))
    
    DF = pd.DataFrame(Prob)
    return DF
    
def ProbDist2(OriAA,Residues,MutAA,Count,df1,RefGene):# Average for each sites
    d = {'Ori':OriAA,'Residues': Residues, \
         'MutAA': MutAA,'Count': Count};  df = pd.DataFrame(d)
    df = pd.merge(df, df1, on ='Residues',how ='left') #df1 contain residues and Ncount column
    
    Count = df["Count"].tolist(); NCount = df["N.Count"].tolist()
    Ori = df["Ori"].tolist(); Mut = df["MutAA"].tolist()
    Residues = df["Residues"].tolist()
    
    Prob,OriMut = [],[]
    FinProb,Residues2 = [],[]
    for i in range(len(Count)):
        Prob.append(round(Count[i]/NCount[i],5))
        OriMut.append(Ori[i]+Mut[i])# OriMut ==> to identify different type permutation
    
    df["OriMut"] = OriMut;df["Prob"] = Prob
    df2 = df.filter(['OriMut','Prob'], axis=1)
    df2 = df2.groupby(['OriMut']).sum()#Numerator
    df2 = df2.reset_index(level=['OriMut']) #To revive the previous summed up OriMut
    
    OriMut = df2["OriMut"].tolist()
    Residues2.extend([i[:1] for i in OriMut])
    df2["Residues"] = Residues2
    
    
    df3 = df2.filter(['Residues','Prob'], axis=1)
    df3 = df3.groupby(['Residues']).sum()#Denumerator
    df3 = df3.reset_index(level=['Residues']) #To revive the previous summed up Residues
    # 
    
    df2 = pd.merge(df2, df3, on ='Residues',how ='left')#final probability
    Num = df2["Prob_x"].tolist();Den = df2["Prob_y"].tolist();Permut = df2["OriMut"].tolist()

    for i in range(len(Num)):
        FinProb.append(round(Num[i]/Den[i],5))
    
    d4 = {'Permut':Permut,'Prob': FinProb};  df4 = pd.DataFrame(d4)
    d5 = {'Permut':RefGene};  df5 = pd.DataFrame(d5)
    df5 = pd.merge(df5, df4, on ='Permut',how ='left')#final data frame
    df5 =df5.fillna(0)
    print(df5)
    
    Prob = df5["Prob"].tolist()
    x = int(np.sqrt(len(RefGene))); y = int(np.sqrt(len(RefGene)))
    Prob = np.array(Prob).reshape((x,y))
    
    DF = pd.DataFrame(Prob)
    return DF      
    
def TailSites(dfRBD,z):
    
    Residue = dfRBD["Sites"].tolist()
    NAA = dfRBD["N.AA"].tolist()
    Thres = math.ceil (st.mean(NAA) + st.stdev(NAA)) #Determining the Threshold
    
    NAASel = []
    for i in range(len(NAA)):
        if int(NAA[i]) >= Thres:
            NAASel.append(int(NAA[i]))
        else:
            NAASel.append(0)
    dSitesRBD= {'Sites': Residue,'N.AA': NAASel}; dfSitesRBD = pd.DataFrame(dSitesRBD) #df for the selected sites, 0 for non selected sites
    # dfSitesRBD["space"] = ""
    
    dfRBD = dfRBD.loc[dfRBD["N.AA"] >= Thres]#select the sites with NAA larger than threshold
    SitesRBD = dfRBD["Sites"].tolist()
    # print("N.Sites =", len(SitesRBD), "Threshold =", Thres)
    
    NAA, NAAfreq = count_dups(NAA)
    dOptionRBD = {'Option': NAA,str(z): NAAfreq}; dfOptionRBD = pd.DataFrame(dOptionRBD)
    dfOptionRBD["Option"] = dfOptionRBD["Option"].astype(int)
    
    dSelectedRBD = {str(z): SitesRBD}; dfSelectedRBD = pd.DataFrame(dSelectedRBD) #df that only contains heavy tail sites
    # dfSelectedRBD["space"] = ""
    
    return SitesRBD, dfSitesRBD, dfSelectedRBD, dfOptionRBD
    
def MutationType(df,ref0,z):
    
    FullResidues = []
    for m in range(len(ref0)-1):
        w = ref0[m]+str(m+1)
        FullResidues.append(str(w))
    
    SNP, DNP, TNP, DEL, Diag, NoDiag = [], [], [], [], [], []
    
    
    Residue1= df["Residues"].tolist()
    NC = df["NC"].tolist()
    Residue1, NC_List = concat_dups(Residue1,NC)#NucleotideChange Counts
    
    Residue2= df["Residues"].tolist()
    MutationType = df["MutationType"].tolist()
    Residue2, MutationType = concat_dups(Residue2,MutationType)#Diag/NoDiag Change
    
    
    d1 = {'Sites': FullResidues}; df1 = pd.DataFrame(d1)
    d2 = {'Sites': Residue1,'NC_List': NC_List}; df2 = pd.DataFrame(d2)#NucleotideChange Counts
    df3 = pd.merge(df1, df2, on ='Sites',how ='left')
    
    d30 = {'Sites': Residue2,'MutationType': MutationType}; df30 = pd.DataFrame(d30)#Diag/NoDiag Change
    df3 = pd.merge(df3, df30, on ='Sites',how ='left')
    
    df3 = df3.fillna(0)
    
    NC_List= df3["NC_List"].tolist() #List of Nucleotide Change .count('a')
    MutationType= df3["MutationType"].tolist() #List of MutationType
    
    for i in range(len(NC_List)):
        snp = str(NC_List[i]).count('1')
        dnp = str(NC_List[i]).count('2')
        tnp = str(NC_List[i]).count('3')
        Del = str(NC_List[i]).count('Del')
        diag = str(MutationType[i]).count('Small')
        nodiag = str(MutationType[i]).count('Big')
        
        SNP.append(snp)
        DNP.append(dnp)
        TNP.append(tnp)
        DEL.append(Del)
        Diag.append(diag)
        NoDiag.append(nodiag)
        
    df3["SNP"] = SNP
    df3["DNP"] = DNP
    df3["TNP"] = TNP
    df3["Del"] = DEL
    df3["Diagonal"] = Diag
    df3["Non-Diagonal"] = NoDiag
    
    dfSNP = df3[['Sites', 'SNP']].copy()
    dfDNP = df3[['Sites', 'DNP']].copy()
    dfTNP = df3[['Sites', 'TNP']].copy()
    dfDel = df3[['Sites', 'Del']].copy()
    dfDiag = df3[['Sites', 'Diagonal']].copy()
    dfNoDiag = df3[['Sites', 'Non-Diagonal']].copy()
        
    return dfSNP, dfDNP, dfTNP, dfDel, dfDiag, dfNoDiag

def Selection(OriAA,Pos,Residues,MutAA,Count,WL_Count,Mutation):
    Col1, Col2, Col3, Col4, Col5, Col6, Col7 = [],[],[],[],[],[],[]
    for i in range(len(Pos)):
        if MutAA[i] != '-':
            # print(Residues[i])
            Col1.append(OriAA[i])
            Col2.append(Pos[i])
            Col3.append(Residues[i])
            Col4.append(MutAA[i])
            Col5.append(Count[i])
            Col6.append(WL_Count[i])
            Col7.append(Mutation[i])

            
    return Col1, Col2, Col3, Col4, Col5, Col6, Col7

###############################################################################################
require_cols1 = [0,1]#Df of AAtype
dfVarRef = pd.read_excel('InputSeq.xlsx', sheet_name ='VarRef', header = 0,usecols = require_cols1)

VarRef = dfVarRef["SeqMSA"].tolist()
VarIndex = 1 


###############################################################################################
file_input = "Reference0Prot.txt"
g = open(file_input,'r')
f0 = g.readlines()
l0 = ref0 = VarRef[0] # Selecting reference seq, 0:Wuhan; 1:Alpha; etc.
length0 = len(ref0)-1

oriAA,oriAARBD,oriAANTD = [],[],[]

###############################################################################################
for m in range(len(ref0)):#List of Originial Residues (1273 AA)
    w = ref0[m]+str(m+1)
    oriAA.append(str(w))

for k in range(325):#List of Originial Residues NTD (1273 AA)
    w = ref0[k]+str(k+1)
    oriAANTD.append(str(w))

for l in range(325,526):#List of Originial Residues RBD (1273 AA)
    w = ref0[l]+str(l+1)
    oriAARBD.append(str(w))

dSites = {'Sites': oriAA}; dfSites = pd.DataFrame(dSites)

#Initial Frame for Basic Data Information
dfPassNAA = dfSites 
dfPassWL_AA = dfSites
dfPassNCount = dfSites
dfPassWL_Count = dfSites

dfPassSNP = dfSites
dfPassDNP = dfSites
dfPassTNP = dfSites
dfPassDel = dfSites
dfPassDiag = dfSites
dfPassNoDiag = dfSites

#Initial Frame for Domain Sites
dSitesAll0 = {'Sites': oriAA}; dfSitesAll0 = pd.DataFrame(dSitesAll0)
dSitesNTD0 = {'Sites': oriAANTD}; dfSitesNTD0 = pd.DataFrame(dSitesNTD0)
dSitesRBD0 = {'Sites': oriAARBD}; dfSitesRBD0 = pd.DataFrame(dSitesRBD0)

#Initial Frame for Option Dist for each Domains
dChannelAll = {'Option': list(range(12))}; dfChannelAll = pd.DataFrame(dChannelAll)
dChannelRBD = {'Option': list(range(12))}; dfChannelRBD = pd.DataFrame(dChannelRBD)
dChannelNTD = {'Option': list(range(12))}; dfChannelNTD = pd.DataFrame(dChannelNTD)
###############################################################################################

loc = ("D:/Research/2 CRISPR/Codes/TO/WL-likePlot/InputSeq.xlsx")#Directory File

wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(0)
l_r = sheet.nrows
l_c = sheet.ncols

###############################################################################################

require_cols = [0,1,4,5]#Df of AAtype
dfTypeAA = pd.read_excel('InputSeq.xlsx', sheet_name ='CodonPermutation', header = 0, \
                   usecols = require_cols)

# TypeAA= dfTypeAA["TypeAA"].tolist()
# TypeAA,TypeAAFreq = count_dups(TypeAA)
# dTypeAA0 = {'TypeAA': TypeAA}; dfTypeAA0 = pd.DataFrame(dTypeAA0) 

###############################################################################################

N = 25 #control the data selection ( How many months that you want to use)
ProtLen = 1273
Numb = range(1,ProtLen+1)
dfPass = pd.DataFrame(Numb) #marker for each row in df when printed to excel
FileProcessing = "SingleSites-TE-0204-noDel-AllUnique-"+str(N) #-NoDel #the ouput file

dDeviation = {'Deviation': list(range(100))}; dfDeviation = pd.DataFrame(dDeviation)
dVariant = {'VarIndex': list(range(12))}; dfVariant = pd.DataFrame(dVariant)

###############################################################################################

for z in range(1,N+1):# z is
    cut = 1 #Cutoff for sequences occurence
    seq,seq0,Variant = [], [], []
    
    # ref = VarRef[1] #Selecting the reference seq
    
    # if z >= 9: 
    #     ref0 = ref
    # else:
    #     ref0 = ref0
    
    for i in range(1,l_r):
        t = int(sheet.cell_value(i,0))
        Mut = sheet.cell_value(i,1)
        FreqRed = sheet.cell_value(i,2)
        Var = sheet.cell_value(i,3)
        
        # if (t != 1 ) and int(FreqRed) >= cut:
        #     if t in range(z-1,z+1):
        #         seq.append(Mut)
        # else:
        #     if t == z and int(FreqRed) >= cut:
        #         seq.append(Mut)
        
        if t == z :
            seq.append(Mut)
            Variant.append(Var)

        # if (t >= z and t < 21) and int(FreqRed) >= cut:
        #     seq.append(Mut)

        # if t == z:
        #     seq.append(Mut)
        
        # if t <= z and int(FreqRed) >= cut:
        #     seq.append(Mut)
    # print('\n')
    
    # seq.append(VarRef[0])
    # seq.append(VarRef[0])
    # Variant.append(0)
    # Variant.append(0)
    
    
    print('iteration -',z, ' Total Seq = ',len(seq), "CutoffSeq =", cut)

    OriAA,Pos,Residues,MutAA,Count,WL_Count,NC = raw_WL(ref0,ProtLen,seq) #PostProcessingDataFrame
    OriAA,Pos,Residues,MutAA,Count,WL_Count,NC = Selection(OriAA,Pos,Residues,MutAA,Count,WL_Count,NC) #Creating 
    df,dfFinal,dfNAA,dfNCount,dfWL_AA,dfWL_Count = WL(OriAA,Pos,Residues,MutAA,Count,WL_Count,NC,ref0,dfTypeAA)#df3 = The list of complete spike
    # dfSNP, dfDNP, dfTNP, dfDel, dfDiag, dfNoDiag = MutationType(df,ref0,z) #determining Mutation Type SNP,DNP,TNP,Diag,NoDiag
    
    ############################################################################################### 
    """ Creating DeviationCounting """
    
    
    DevSeq = []
    for j in range(len(seq)):
        Dev = 0 #counting N.Mut
        for k in range(len(seq[j])):
            if seq[j][k] != ref0[k]: Dev += 1
        DevSeq.append(Dev)
    
    Deviation, DevFreq = count_dups(DevSeq)
    # for i in range(len(DevFreq)):
    #     DevFreq[i] = DevFreq[i]/len(seq)
    dDev = {'Deviation': Deviation,str(z): DevFreq}; dfDev = pd.DataFrame(dDev)
    # print(dfDev)
    dfDeviation = pd.merge(dfDeviation, dfDev, on ='Deviation',how ='outer')
    dfDeviation = dfDeviation.fillna(0)
    
    
    ############################################################################################### 
    """ Creating Variant Counting """
    
    
    VarIndex, VarFreq = count_dups(Variant)
    dVar = {'VarIndex': VarIndex,str(z): VarFreq}; dfVar = pd.DataFrame(dVar)
    dfVariant = pd.merge(dfVariant, dfVar, on ='VarIndex',how ='outer')
    dfVariant = dfVariant.fillna(0)
    ############################################################################################### 
    """ Creating dfNAA, dfWL_AA, dfNCount, dfWL_Count """
    
    dfNAA = dfNAA.rename(columns={'N.AA':str(z)})
    dfPassNAA = pd.merge(dfPassNAA,dfNAA, on ='Sites',how ='left')
    
    dfWL_AA = dfWL_AA.rename(columns={'WL-AA':str(z)})
    dfPassWL_AA = pd.merge(dfPassWL_AA,dfWL_AA, on ='Sites',how ='left')
    
    dfNCount = dfNCount.rename(columns={'N.Count':str(z)})
    dfPassNCount = pd.merge(dfPassNCount,dfNCount, on ='Sites',how ='left')
    
    dfWL_Count  = dfWL_Count.rename(columns={'WL-Count':str(z)})
    dfPassWL_Count = pd.merge(dfPassWL_Count,dfWL_Count, on ='Sites',how ='left')

    ###############################################################################################
    """ Selecting All Sites """
    
    dfAll = dfFinal #Selecting All Sites
    SitesAll, dfSitesAll ,dfSelectedAll, dfOptionAll = TailSites(dfAll,z)
    
    dfChannelAll = pd.merge(dfChannelAll, dfOptionAll, on ='Option',how ='outer')
    dfChannelDistAll = dfChannelAll.fillna(0)
    dfChannelDistAll["space"] = ""
    
    dfSitesAll0 = pd.merge(dfSitesAll0, dfSitesAll, on ='Sites',how ='left') #df for the selected sites, 0 for non selected sites
    dfSitesAll0 = dfSitesAll0.rename(columns={'N.AA':str(z)})
    dfSitesAll0 = dfSitesAll0.fillna(0)
    
    ###############################################################################################
    """ Selecting NTD Sites """
    
    dfNTD = dfFinal.loc[dfFinal["Pos"] <= 324] #Selecting NTD Sites
    SitesNTD, dfSitesNTD ,dfSelectedNTD, dfOptionNTD = TailSites(dfNTD,z)
    
    dfChannelNTD = pd.merge(dfChannelNTD, dfOptionNTD, on ='Option',how ='outer')
    dfChannelDistNTD = dfChannelNTD.fillna(0)
    dfChannelDistNTD["space"] = ""
    
    dfSitesNTD0 = pd.merge(dfSitesNTD0, dfSitesNTD, on ='Sites',how ='left') #df for the selected sites, 0 for non selected sites
    dfSitesNTD0 = dfSitesNTD0.rename(columns={'N.AA':str(z)})
    dfSitesNTD0 = dfSitesNTD0.fillna(0)

    ###############################################################################################
    """ Selecting RBD Sites """
    
    dfRBD = dfFinal[(dfFinal['Pos'] >= 325) & (dfFinal['Pos'] <= 525)]
    SitesRBD, dfSitesRBD ,dfSelectedRBD, dfOptionRBD = TailSites(dfRBD,z)
    
    dfChannelRBD = pd.merge(dfChannelRBD, dfOptionRBD, on ='Option',how ='outer')
    dfChannelDistRBD = dfChannelRBD.fillna(0)
    dfChannelDistRBD["space"] = ""
    
    dfSitesRBD0 = pd.merge(dfSitesRBD0, dfSitesRBD, on ='Sites',how ='left')
    dfSitesRBD0 = dfSitesRBD0.rename(columns={'N.AA':str(z)})
    dfSitesRBD0 = dfSitesRBD0.fillna(0)

    ###############################################################################################
    """ Creating SNP,DNP,TNP,Del,Diag,No-Diag data frame """
    
    # dfSNP = dfSNP.rename(columns={'SNP':str(z)})
    # dfPassSNP = pd.merge(dfPassSNP,dfSNP, on ='Sites',how ='left') #dfSNP, dfDNP, dfTNP, dfDel, dfDiag, dfNoDiag
    
    # dfDNP = dfDNP.rename(columns={'DNP':str(z)})
    # dfPassDNP = pd.merge(dfPassDNP, dfDNP, on ='Sites',how ='left')
    
    # dfTNP = dfTNP.rename(columns={'TNP':str(z)})
    # dfPassTNP = pd.merge(dfPassTNP, dfTNP, on ='Sites',how ='left')
    
    # dfDel  = dfDel.rename(columns={'Del':str(z)})
    # dfPassDel = pd.merge(dfPassDel, dfDel, on ='Sites',how ='left')
    
    # dfDiag  = dfDiag.rename(columns={'Diagonal':str(z)})
    # dfPassDiag = pd.merge(dfPassDiag, dfDiag, on ='Sites',how ='left')
    
    # dfNoDiag  = dfNoDiag.rename(columns={'Non-Diagonal':str(z)})
    # dfPassNoDiag = pd.merge(dfPassNoDiag, dfNoDiag, on ='Sites',how ='left')
    

###############################################################################################
""" Print Out Basic Information: dfNAA, dfWL_AA, dfNCount, dfWL_Count """

dfPassNAA["space"] = ""
dfPassWL_AA["space"] = ""
dfPassNCount["space"] = ""
dfPassWL_Count["space"] = ""
df1 = pd.concat([dfPassNAA,dfPassWL_AA,dfPassNCount,dfPassWL_Count], axis =1)
file_name1 = FileProcessing+'.xlsx'#
df1.to_excel(file_name1, sheet_name='BasicData', index=None)


###############################################################################################

path = r"D:/Research/2 CRISPR/Codes/TO/WL-likePlot/" + file_name1
book = load_workbook(path)
writer = pd.ExcelWriter(path, engine = 'openpyxl')
writer.book = book

###############################################################################################

""" Print Out Deviation and Variant Counting"""
dfDeviation["space"] = ""
dfDeviation = pd.concat([dfDeviation,dfVariant], axis =1)
dfDeviation.to_excel(writer, sheet_name='Dev-Var',index=None)

###############################################################################################

""" Print Out AllSites"""
dfChannelDistAll = dfChannelDistAll.rename(columns={'Option':'Option-All'})
dfSitesAll0["space"] = ""
dfSitesFinalAll = dfSitesAll0.rename(columns={'Sites':'All-Sites'})
df2 = pd.concat([dfChannelDistAll,dfSitesFinalAll], axis =1)
df2.to_excel(writer, sheet_name='All',index=None)

###############################################################################################

""" Print Out NTD"""
dfChannelDistNTD = dfChannelDistNTD.rename(columns={'Option':'Option-NTD'})
dfSitesNTD0["space"] = ""
dfSitesFinalNTD = dfSitesNTD0.rename(columns={'Sites':'NTD-Sites'})
df2 = pd.concat([dfChannelDistNTD,dfSitesFinalNTD], axis =1)
df2.to_excel(writer, sheet_name='NTD',index=None)

###############################################################################################

""" Print Out RBD"""
dfChannelDistRBD = dfChannelDistRBD.rename(columns={'Option':'Option-RBD'})
dfSitesRBD0["space"] = ""
dfSitesFinalRBD = dfSitesRBD0.rename(columns={'Sites':'RBD-Sites'})
df3 = pd.concat([dfChannelDistRBD,dfSitesFinalRBD], axis =1)
df3.to_excel(writer, sheet_name='RBD',index=None)

###############################################################################################
""" Print Out PostProcessing"""
df = pd.concat([df,dfFinal], axis =1)
df.to_excel(writer, sheet_name='PostProcessing',index=None)#writer

###############################################################################################
""" Print Out MutationType"""
dfPassSNP["space"] = ""
dfPassSNP = dfPassSNP.rename(columns={'Sites':'SNP'})
dfPassSNP.to_excel(writer, sheet_name='SNP',index=None)

dfPassDNP["space"] = ""
dfPassDNP = dfPassDNP.rename(columns={'Sites':'DNP'})
dfPassDNP.to_excel(writer, sheet_name='DNP',index=None)

dfPassTNP["space"] = ""
dfPassTNP = dfPassTNP.rename(columns={'Sites':'TNP'})
dfPassTNP.to_excel(writer, sheet_name='TNP',index=None)

dfPassDel["space"] = ""
dfPassDel = dfPassDel.rename(columns={'Sites':'Del'})
dfPassDel.to_excel(writer, sheet_name='Deletion',index=None)

dfPassDiag["space"] = ""
dfPassDiag = dfPassDiag.rename(columns={'Sites':'Diagonal'})
dfPassDiag.to_excel(writer, sheet_name='Diagonal',index=None)

dfPassNoDiag["space"] = ""
dfPassNoDiag = dfPassNoDiag.rename(columns={'Sites':'Non-Diagonal'})
dfPassNoDiag.to_excel(writer, sheet_name='Non-Diagonal',index=None)

dfMutType = pd.concat([dfPassSNP, dfPassDNP, dfPassTNP, dfPassDel, dfPassDiag, dfPassNoDiag], axis =1)
dfMutType.to_excel(writer, sheet_name='MutType',index=None)

##############################################################################################

writer.save()

print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')
