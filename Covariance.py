# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 09:30:46 2021

@author: Hasan
"""

import xlrd
import pandas as pd
import time
import numpy as np
import statistics as st
import math
from datetime import datetime
from scipy.linalg import svdvals

def insert_n(f1):
    f2=[]
    for i in range(len(f1)):
        f2.append(f1[i])
        f2.append("\n")
    return f2

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

def Covariance(oriAA,seq):#,dfAAPos1ref,dfAAPos2ref,dfPairAAref ProtLen,
            
    ResiduePair, CovScore = [],[]
    
    for i in range(ProtLen):#Loop every position
        for j in range(ProtLen):
            if i !=j  and j > i:
                AAPos1, AAPos2, PairAA = [],[],[]  #to store AA for each seq
                Cov12 = []
                
                Pos1 = i+1
                Pos2 = j+1
                
                
                # print("Start =", current_time)
                
                for l in range(len(seq)):#Loop every seq
                    seqj = seq[l] #assign a seq
                    AAPos1.append(seqj[Pos1-1]) #record the AA of the seq on the 1st position
                    AAPos2.append(seqj[Pos2-1]) #record the AA of the seq on the 2nd position
                    PairAA.append(seqj[Pos1-1]+seqj[Pos2-1]) #record the pair AA of the seq on the 2nd position
                
                
                AAPos1,AAPos1Freq = count_dups(AAPos1)
                dAAPos1 = {'AA1': AAPos1, 'FreqAA1': AAPos1Freq}; dfAAPos1  = pd.DataFrame(dAAPos1)#Pos1
                dfAAPos1Final = pd.merge(dfAAPos1ref, dfAAPos1, on ='AA1',how ='left')
                dfAAPos1Final = dfAAPos1Final.fillna(0)
                FreqAA1 = dfAAPos1Final["FreqAA1"].tolist()
                
                AAPos2,AAPos2Freq = count_dups(AAPos2)
                dAAPos2 = {'AA2': AAPos2, 'FreqAA2': AAPos2Freq}; dfAAPos2  = pd.DataFrame(dAAPos2)#Pos2
                dfAAPos2Final = pd.merge(dfAAPos2ref, dfAAPos2, on ='AA2',how ='left')
                dfAAPos2Final = dfAAPos2Final.fillna(0)
                FreqAA2 = dfAAPos2Final["FreqAA2"].tolist()
                
                PairAA,PairAAFreq = count_dups(PairAA)
                dPairAA = {'PairAA': PairAA, 'FreqPairAA': PairAAFreq}; dfPairAA  = pd.DataFrame(dPairAA)#PairAA
                dfPairAAFinal = pd.merge(dfPairAAref, dfPairAA, on ='PairAA',how ='left')
                dfPairAAFinal = dfPairAAFinal.fillna(0)
                FreqPairAA = dfPairAAFinal["FreqPairAA"].tolist()
                
                for k in range(len(FreqPairAA)):
                    c = ( FreqPairAA[k]/len(seq) ) - ( (FreqAA1[k] * FreqAA2[k]) / len(seq)**2 )
                    Cov12.append( c )
                        
                Cov12 = np.asarray(Cov12)
                Cov12 = Cov12.reshape(21, 21)
                
                w = svdvals(Cov12)#applying SVD
                p = 0
                
                for l in range(0,len(w)):
                    p += (w[l])**2
                
                p = np.sqrt(p)
                
                ResiduePair.append(oriAA[i]+'-'+oriAA[j])
                CovScore.append(p)
                
                # print(oriAA[i],'-',oriAA[j],'=',p)
    
    dCovScore = {'Pair': ResiduePair, 'CovScore': CovScore}; dfCovScore  = pd.DataFrame(dCovScore)
    
    return dfCovScore 

###############################################################################################
file_input = "Reference0Prot.txt"
g = open(file_input,'r')
f0 = g.readlines()
l0 = ref0 = f0[0]
length0 = len(ref0)-1

oriAA = []
for m in range(len(ref0)-1):#List of Originial Residues (1273 AA)
    w = ref0[m]+str(m+1)
    oriAA.append(str(w))

###############################################################################################
loc = ("D:/Research/2 CRISPR/Codes/TO/WL-likePlot/InputSeq.xlsx")#Directory File
wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(0)
l_r = sheet.nrows
l_c = sheet.ncols

###############################################################################################

ProtLen = 1273
list_AA0 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']
list_AA, list_AA1, list_AA2 = [], [], []

for aa in list_AA0: #defining the pair of AA
    for bb in list_AA0:
        list_AA.append(aa+bb)
        list_AA1.append(aa)
        list_AA2.append(bb)
        
dAAPos1ref = {'AA1': list_AA1}; dfAAPos1ref  = pd.DataFrame(dAAPos1ref)
dAAPos2ref = {'AA2': list_AA2}; dfAAPos2ref  = pd.DataFrame(dAAPos2ref)
dPairAAref = {'PairAA': list_AA}; dfPairAAref  = pd.DataFrame(dPairAAref)


###############################################################################################
N = 21 #control the data selection

start_time = time.process_time()
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Start =", current_time)

for z in range(14,N+1):# z is
    cut = 1 #Cutoff for sequences occurence
    seq = []
    for i in range(1,l_r):
        t = int(sheet.cell_value(i,0))
        Mut = sheet.cell_value(i,1)
        FreqRed = sheet.cell_value(i,2)
                
        if t == z and int(FreqRed) >= cut:
            seq.append(Mut)

    print('\niteration -Month',z, ' Total Seq = ',len(seq), "CutoffSeq =", cut)
    dfCovScore = Covariance(oriAA,seq)    #,dfAAPos1ref,dfAAPos2ref,dfPairAAref
    dfCovScore.to_excel('Covariance-Month'+str(z)+'-Cutoff'+str(cut)+'.xlsx',index=None)#writer
    

    print('TIME TAKEN: ' + str(time.process_time() - start_time) + 's')
    
    
"""
    # seq = insert_n(seq)
    # with open('SeqInputCoV-Month'+str(z)+'.txt', 'w') as f:
    #     for item in seq:
    #         f.write("%s" % item)

"""   
        












