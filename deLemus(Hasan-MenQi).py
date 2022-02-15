# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 13:57:30 2022

@author: Hasan (In Colaboration with MenQi)
"""

import xlrd
import pandas as pd
import time
import numpy as np
import statistics as st
import math
from datetime import datetime
from scipy.linalg import svd
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

def split(word):
    return [char for char in word] #split a string into each word

def ScoringFunc(ref0,seq,dfRef1):
    Mut = []
    seq = split(seq)
    ref0 = split(ref0)
    
    for m in range(len(ref0)):
        w = ref0[m]+seq[m]
        Mut.append(str(w))
    
    dseq = {'XY': Mut}; dfseq = pd.DataFrame(dseq)
    dfseq =  pd.merge(dfseq, dfRef1, on ='XY',how ='left')
    Score = dfseq["ScoreBinary"].tolist()
    
    return Score

def SVDfunc(A,kk):
    U, s, VT = svd(A)
    r = A.shape[0]
    rs = s.shape[0]
    Au=0.0*A
    
    Sum_s = 0.0
    for i in range(len(s)): Sum_s += s[i]**2
    for j in range(len(s)): s[j] = s[j]**2
    s = s[:6]
    # print(s)
    s = s/Sum_s
    # print(Sum_s)
    # print(s)
    
    for b in range(rs):
        Auu = 0.0
        for k in range(r):
            Auu = Auu + U[k,b]*A[k,:]
        Au[b,:] = Auu
    
    Au = Au[:kk,:]
    r,c = Au.shape
    
    SVDScore = []
    for y in range(c):#loop over columns
        c = 0.0
        for x in range(r):#loop over rows
            c += (Au[x,y]**2)
        SVDScore.append(c**0.5)
    
    return SVDScore,s

def LogScale(df):
    df = df.astype(float)
    y = df.to_numpy()
    for i in range(len(y)):
        for j in range(len(y[0])):
            y0 = y[i,j]
            y[i,j] = math.log10(float(y0) + 1.0)
            # print(y0,y[i,j])
    df = pd.DataFrame(y)
    return df

###############################################################################################
loc = ("D:/Research/2 CRISPR/Codes/TO/WL-likePlot/InputSeq.xlsx")#Directory File

wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(0)
l_r = sheet.nrows
l_c = sheet.ncols

###############################################################################################
require_cols1 = [0,1]#Df of AAtype
dfVarRef = pd.read_excel('InputSeq.xlsx', sheet_name ='VarRef', header = 0, \
                    usecols = require_cols1)

VarRef = dfVarRef["SeqMSA"].tolist()
VarIndex = 1 
ref0 = VarRef[0]

oriAA = []
for m in range(len(ref0)):#List of Originial Residues (1273 AA)
    w = ref0[m]+str(m+1)
    oriAA.append(str(w))

###############################################################################################
require_cols = [0,1,4,5]#Df of AAtype
dfMutTypeAA = pd.read_excel('InputSeq.xlsx', sheet_name ='CodonPermutation', header = 0, \
                   usecols = require_cols)

dDeviation = {'Deviation': list(range(100))}; dfDeviation = pd.DataFrame(dDeviation)
dVariant = {'VarIndex': list(range(12))}; dfVariant = pd.DataFrame(dVariant)
dSites = {'Sites': oriAA}; dfPassNAA = pd.DataFrame(dSites)
dChannelAll = {'Option': list(range(12))}; dfChannelAll = pd.DataFrame(dChannelAll)
dDev = {'Deviation': list(range(100))}; dfMeanDev = pd.DataFrame(dDev)
dSAP= {'SAP': list(range(12))}; dfMeanSAP = pd.DataFrame(dSAP)


require_cols = [2,3,4,5]##define the require columns
dfTypeAA = pd.read_excel('InputSeq.xlsx', sheet_name ='TypeAAPermutation', header = 0, \
                   usecols = require_cols)
dfRef1 = dfTypeAA[['XY', 'ScoreBinary']]#Binary
dSites = {'Sites': range(1,1273+1)}; dfSites = pd.DataFrame(dSites)

###############################################################################################
N = 25 #control the data selection
kk = 4 #control the number of eigenvalue
ProtLen = 1273
MeanDev,MeanSAP = [],[]
Eigen = []
iter_range = range(1,N+1)
for z in iter_range:# z is
    seq,seq0,Variant = [], [], []
    
    for i in range(1,l_r):
        t = int(sheet.cell_value(i,0))
        Mut = sheet.cell_value(i,1)
        FreqRed = sheet.cell_value(i,2)
        Var = sheet.cell_value(i,3)
          
        if t == z :  #and int(FreqRed) >= cut and Var >= 0
            seq.append(Mut)
            Variant.append(Var)
    # require_cols = [0,1,2,3]#Df of AAtype
    # dfInput = pd.read_excel('InputSeq.xlsx', sheet_name ='MonthlyUnique1812', header = 0, \
    #                usecols = require_cols)
    # dfInput = dfInput.loc[dfInput["MonthIndex"] == z]
    # seq = dfInput["SeqMSA"].tolist()
    # Variant = dfInput["Detect VoC"].tolist()
    
    print('iteration -',z, ' Total Seq = ',len(seq))
    OriAA,Pos,Residues,MutAA,Count,WL_Count,NC = raw_WL(ref0,ProtLen,seq)
    OriAA,Pos,Residues,MutAA,Count,WL_Count,NC = Selection(OriAA,Pos,Residues,MutAA,Count,WL_Count,NC)
    df,dfFinal,dfNAA,dfNCount,dfWL_AA,dfWL_Count = WL(OriAA,Pos,Residues,MutAA,Count,WL_Count,NC,ref0,dfMutTypeAA)
    
    ############################################################################################### 
    """ Creating DeviationCounting """
    DevSeq = []
    for j in range(len(seq)):
        Dev = 0 #counting N.Mut
        for k in range(len(seq[j])):
            if seq[j][k] != ref0[k]: Dev += 1
        DevSeq.append(Dev)
    
    Deviation, DevFreq = count_dups(DevSeq)
    for i in range(len(DevFreq)):
        DevFreq[i] = DevFreq[i]/len(seq)
    dDev = {'Deviation': Deviation,str(z): DevFreq}; dfDev = pd.DataFrame(dDev)
    dfDeviation = pd.merge(dfDeviation, dfDev, on ='Deviation',how ='outer')
    dfDeviation = dfDeviation.fillna(0)
    
    dfMeanDev[str(z)]= dfDeviation["Deviation"]*dfDeviation[str(z)]
    MeanDev.append( dfMeanDev[str(z)].sum() )
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
    
    ###############################################################################################
    """ Selecting All Sites """
    dfAll = dfFinal #Selecting All Sites
    SitesAll, dfSitesAll ,dfSelectedAll, dfOptionAll = TailSites(dfAll,z)
    dfChannelAll = pd.merge(dfChannelAll, dfOptionAll, on ='Option',how ='outer')
    dfChannelDistAll = dfChannelAll.fillna(0)
    
    dfMeanSAP[str(z)]= dfChannelDistAll["Option"]*dfChannelDistAll[str(z)]
    MeanSAP.append( dfMeanSAP[str(z)].sum()/ProtLen )
    
    ###############################################################################################
    """ CreatingInputSVD MAtrix and SVD Calculation """
    SeqArray = [] #SeqBinary = Convert each seq to be binary mutation landscape, SeqArray = the total binary landscape of selected time range for all seq
    for j in range(len(seq)):#loop over sequences
    
        Seq0 = seq[j]#assign every seq
        w = ScoringFunc(ref0,Seq0,dfRef1)
        SeqArray.append(w)
    
    SeqArray = np.array(SeqArray)
    A = SeqArray.astype(float)
    A = np.nan_to_num(A)
    SVDScore,s = SVDfunc(A,kk)
    Eigen.append(s)
    dSVD = {'Sites': range(1,1273+1),str(z): SVDScore}; dfSVD = pd.DataFrame(dSVD)
    dfSites = pd.merge(dfSites,dfSVD, on ='Sites',how ='left')
    
    
    print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')

dfSAP = dfPassNAA.iloc[:, 1:len(dfPassNAA.columns) ]# [range(1,len(dfPassNAA.columns))]
dfSVD = dfSites.iloc[:, 1:len(dfSites.columns) ]
ddeLemus = {'Pos': range(1,ProtLen+1)}; dfdeLemus = pd.DataFrame(ddeLemus)
dfdeLemus = pd.concat([dfdeLemus,dfSAP.mul(dfSVD)], axis =1)
dfdeLemusArray = LogScale(dfdeLemus)
deLemus = dfdeLemusArray.to_numpy()

df0 = pd.DataFrame(deLemus)
file_name = 'deLemus0204-AllUnique'+str(N)+'-VarAdjustedTrial.xlsx'#-Cumulative-2LastMonth+str(cut)-Deletion
df0.to_excel(file_name, sheet_name='deLemus')
###############################################################################################
path = r"D:/Research/2 CRISPR/Codes/TO/WL-likePlot/" + file_name
book = load_workbook(path)
writer = pd.ExcelWriter(path, engine = 'openpyxl')
writer.book = book


""" Print Out Mean SAP and Deviation"""
dMean = {'Month': range(1,N+1)}; dfMean= pd.DataFrame(dMean)
dfMean['MeanDev'] = MeanDev
dfMean['MeanSAP'] = MeanSAP
dfMean.to_excel(writer, sheet_name='MeanSAP-Dev',index=None)

###############################################################################################
dDynNTDSele = {'Pos': range(1,45+1)}; dfDynNTDSele= pd.DataFrame(dDynNTDSele)
dDynRBDSele = {'Pos': range(1,20+1)}; dfDynRBDSele= pd.DataFrame(dDynRBDSele)
DynNTDSele ,DynRBDSele,NTD,RBD = [], [],[], []
for z in iter_range:
    dfNTD = dfdeLemus[(dfdeLemus['Pos'] < 325) ]
    dfNTD.sort_values(by=str(z), inplace=True,ascending=False )# time sorted data
    dfNTDSelected = dfNTD.head(45)
    NTDSelected = sorted( dfNTDSelected['Pos'].tolist() )
    dfDynNTDSele[str(z)] = NTDSelected
    NTD.extend(NTDSelected)
    
    dfRBD = dfdeLemus[(dfdeLemus['Pos'] >= 325) & (dfdeLemus['Pos'] <= 525)]
    dfRBD.sort_values(by=str(z), inplace=True,ascending=False )# time sorted data
    dfRBDSelected = dfRBD.head(20)
    RBDSelected = sorted( dfRBDSelected['Pos'].tolist() )
    dfDynRBDSele[str(z)] = RBDSelected
    RBD.extend(RBDSelected)

NTD, NTDFreq = count_dups(NTD)
dNTD = {'Sites': NTD,'Freq': NTDFreq}; dfNTD = pd.DataFrame(dNTD)
RBD, RBDFreq = count_dups(RBD)
dRBD = {'Sites': RBD,'Freq': RBDFreq}; dfRBD = pd.DataFrame(dRBD)
""" Print Out SelectedSites"""
dfDynSele = pd.concat([dfDynNTDSele,dfDynRBDSele,dfNTD,dfRBD], axis =1)
dfDynSele.to_excel(writer, sheet_name='SelectedSites',index=None)

###############################################################################################
Eigen = np.array(Eigen)
dfEigen = pd.DataFrame(Eigen)
dfEigen.to_excel(writer, sheet_name='Top6Eigen',index=None)
###############################################################################################
writer.save()

print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')