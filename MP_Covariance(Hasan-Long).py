# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 21:30:22 2021

@author: Hasan
"""

from multiprocessing import Manager
from multiprocessing import Pool
from functools import partial
import numpy as np
import pandas as pd
from scipy.linalg import svdvals
import xlrd

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

def Covariance(ProtLen,oriAA,seq,dfAAPos1ref,dfAAPos2ref,dfPairAAref,output_list,ID):#
            
    # ResiduePair, CovScore = [],[]
    pfinal = 0
    pairRes = 0
    j = ID
    for i in range(ProtLen):#Loop every position
        for j in range(ProtLen):
            if i !=j  and j > i:
                AAPos1, AAPos2, PairAA = [],[],[]  #to store AA for each seq
                Cov12 = []
                
                Pos1 = i+1
                Pos2 = j+1
                
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
                                
                pfinal = np.sqrt(p)
                pairRes = oriAA[i]+'-'+oriAA[j]
                output_list.append([j,pairRes,pfinal])
                

if __name__ == '__main__':
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
    loc = ("/home/mhasana/InputSeqCovariance.xlsx")#Directory File 
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
    
    for z in range(8,N+1):# z is
        cut = 1 #Cutoff for sequences occurence
        seq = []
        for i in range(1,l_r):
            t = int(sheet.cell_value(i,0))
            Mut = sheet.cell_value(i,1)
            FreqRed = sheet.cell_value(i,2)
                    
            
            if t == z and int(FreqRed) >= cut:
                seq.append(Mut)
    
        print('iteration -Month',z, ' Total Seq = ',len(seq), "CutoffSeq =", cut)
        # dfCovScore = Covariance(ProtLen,oriAA,seq,dfAAPos1ref,dfAAPos2ref,dfPairAAref) 

        output_list = []
        output_list = Manager().list() # create shared list
        queue_IDs=range(ProtLen) # get queue ID list ===>should be ProtLen
        print(queue_IDs)
        pool = Pool(processes=160) # input number of core
        func = partial(Covariance,ProtLen,oriAA,seq,dfAAPos1ref,dfAAPos2ref,dfPairAAref,output_list) # change the f1 and f2
        #                ^          ^           ^
        #                |          |           |
        #            function   shared        data
        #            name       list          file
        pool.map(func, queue_IDs) # start multi_processig
        pool.close() # finish multiprocessing
        pool.join() # join the result
        out_df=pd.DataFrame.from_records(output_list,columns=['ID','PairResidue','pfinal'],index='ID')
        out_df.sort_index() # short queue ID
        out_df.to_excel('Covariance-Month'+str(z)+'-Cutoff'+str(cut)+'.xlsx',index=None) # output