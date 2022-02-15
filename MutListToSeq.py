# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 12:26:58 2022

@author: Hasan
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

def strip_n(f1):
    f2=[]
    for i in range(len(f1)):
        f1[i]=f1[i].rstrip('\r\n')
        f2.append(f1[i])
    return f2

def split(word):
    return [char for char in word] #split a string into each word

def OnlyNumb(string):#retrieve number
    word = [char for char in string] #split the string to each character
    w = [s for s in word if s.isdigit()] #keep the number
    new = ""
    for x in w:
        new += x
    
    return int(new)

file_input = "Reference0Prot.txt"
g = open(file_input,'r')
f0 = g.readlines()
f0 = strip_n(f0) #removing '\n'
f0 = split( f0[0] ) #splitting the long string into list of character

require_cols1 = [2,3]
dfMutList = pd.read_excel('COV_Spikes.immediate.2.3(MasterShi-Hong).xlsx', sheet_name ='Sheet2', header = 0,usecols = require_cols1)
MutList = dfMutList['mutation info'].tolist()
ID = dfMutList['ID']
SeqMSA = []

for i in range(len(MutList)):
    
    print(i+1,MutList[i])
    OriSeq = list(f0)
    w = MutList[i].split(";")
    
    for j in range(len(w)):
        NewAA = w[j][-1:]
        Pos = OnlyNumb(w[j]) - 1
        OriSeq[Pos] = NewAA
    
    NewSeq = ''.join(OriSeq)
    SeqMSA.append(NewSeq)

dSeqMSA = {'SeqMSA': SeqMSA}; dfSeqMSA = pd.DataFrame(dSeqMSA)
dfSeqMSA = dfSeqMSA. join(ID)
print(dfSeqMSA)
file_name = 'SeqMSA'+str(len(MutList))+'-178k.xlsx'#-Cumulative-2LastMonth+str(cut)-Deletion
dfSeqMSA.to_excel(file_name, sheet_name='deLemus')


print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')