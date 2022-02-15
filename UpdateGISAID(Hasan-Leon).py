# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 13:01:59 2021

@author: Hasan (Modified by Leon)
"""

import pandas as pd
import numpy as np
import time

from datetime import datetime

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

def insert_n(f1):
    f2=[]
    for i in range(len(f1)):
        f2.append(f1[i])
        f2.append("\n")
    return f2

def split(word):
    return [char for char in word] #split a string into each word

def concat_list(list):
    result= ''
    for element in list:
        result += str(element)
    return result

def count_dups(nums): 
    element = []
    freque = []
    nums = sorted(nums)
    if not nums:
        return element
    running_count = 1
    for i in range(len(nums)-1):
        if nums[i] == nums[i+1]: #
            running_count += 1
        else:
            freque.append(running_count)
            element.append(nums[i])
            running_count = 1
    freque.append(running_count)
    element.append(nums[i+1])
    return element,freque

def concat_dups(criteria,ele): # ele is Amino Acid
    element = []
    freque = []
    if not criteria:
        return element
    
    running_count = ele[0]
    for i in range(len(criteria)-1):
        if criteria[i] == criteria[i+1]:
            running_count += ","+ele[i+1]
        else:
            freque.append(running_count)
            element.append(criteria[i])
            running_count = ele[i+1]

    freque.append(running_count)
    element.append(criteria[i+1])
    return element,freque

def SplitChar(df):
    twodim_list=[]
    for index,line in df.iterrows():
        character_in_line=[]
        line=line.values[0]
        #print(line)
        for character in line:        
            character_in_line.append(character)
        twodim_list.append(character_in_line)
    
    df=pd.DataFrame(twodim_list)
    return df

def countlen(df):
    file_input = df
    g = open(file_input,'r')
    f0 = g.readlines()
    f1 = strip_n(f0)
    result = [] 
    for i in range(len(f1)):
        result.append(len(f1[i]))
        result.append('\n')
    with open('lenseq.txt','w') as f1:#rename this file ----------------------------------------------------------
        for item in result:
            f1.write("%s" % item)

def sumsimilar(criteria,ele): #ele is number
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

def Split(text): #Separate "|", return list
    file_input = text
    g = open(file_input,'r')
    f0 = g.readlines()
    f1 = strip_n(f0)
    result = []
    for i in range(len(f1)):
        f2 = f1[i].split("|")
        result.append(f2[0])
        result.append('\n')
        result.append(f2[1])
        result.append('\n')
    # with open('ProcessedNewUnique-1218.fasta','w') as f5:#rename this file ----------------------------------------------------------
    #     for item in result:
    #         f5.write("%s" % item)
    return result

def SearchDel(Seq):   #Just finding seq with deletion
    g = open(Seq,'r')
    f1 = g.readlines()
    f1 = strip_n(f1)
    dellist = []
    for i in range(len(f1)):
        for j in range(len(f1[i])):
            if f1[i][j] == "-":
                dellist.append(f1[i])
    print("Amount of New Unique Seq With Deletion:", len(dellist))
    with open('Del-Seq1218_100.txt', 'w') as f:#change the number based on the txt file
    	for item in dellist:
    		f.write("%s" % item)

def SearchMoreThan2Del(Seq): #Finding seq with more than 2 deletion
    g = open(Seq,'r')
    f1 = g.readlines()
    f1 = strip_n(f1)
    morethan2 = []
    for i in range(len(f1)):
        dashcount = f1[i].count("-")
        if dashcount >= 3:
            morethan2.append(f1[i])
    print("Amount of Seq w/ >3 Del:", len(morethan2))
    with open('MoreThan2Del-Seq1218_100.txt', 'w') as f:#change the number based on the txt file
     	for item in morethan2:
             f.write("%s" % item)

def Separate(fasta,seq_amount): 
    num = seq_amount*2 ################# times 2 because there are ID and Seq
    file_input = fasta
    g = open(file_input,'r')
    f0 = g.readlines()
    f1 = strip_n(f0)
    j = 0
    while j < len(f1)-num :
        separated_list = []
        for i in range(j,j+num):
            separated_list.append(f1[i])
            separated_list.append('\n')
        filename = fasta[:-6] + "-" + str(seq_amount)+'_'+ str(j) +'.txt'
        with open(filename,'w') as f3:#rename this file ----------------------------------------------------------
            for item in separated_list:
                f3.write("%s" % item)
        print(str(j)+"-"+str(j+num))
        j += num
    length = len(f1)
    print(str(j)+"-"+str(length))
    newlist = []
    for i in range(j,length):
        newlist.append(f1[i])
        newlist.append('\n')
    filename = fasta[:-6] + "-" + str(seq_amount)+'_'+ str(j) +'.txt'
    with open(filename,'w') as f4:#rename this file ----------------------------------------------------------
        for item in newlist:
            f4.write("%s" % item)  
              
def ReadAlign(f1):
    
    length = len(f1) #number of rows in fasta file
    SeqPass, Seq,ID= [],[],[]
    ID.append(f1[0].rstrip('\r\n'))
    b = 0 #initial row which contain ">"
    
    for i in range(1,length) :
        l = f1[i]
        if l[0] != '>':
            SeqPass.append(f1[i].rstrip('\r\n'))
            q = i
        if l[0] == '>': #if f1[i] contain ">"
            ID.append(f1[i].rstrip('\r\n'))
            q = i
            Seq.append(concat_list(SeqPass))
            b = q
            SeqPass = []
    
    SeqPass = strip_n(f1[b+1:length]) #add the last element
    Seq.append(concat_list(SeqPass))
    
    return ID, Seq

def ReadOutFile(outfile):
    file_input = outfile
    g = open(file_input,'r')
    f0 = g.readlines()
    f1 = strip_n(f0)
    result = []
    for i in range(len(f1)):
        result.append(f1[i])
        result.append('\n')
    name = 'OutFileReading' + outfile + '.txt'
    with open(name,'w') as f5:#rename this file ----------------------------------------------------------
        for item in result:
            f5.write("%s" % item)

def PreProcessing(fnew):#remove X AA and unspecific ColDate(month-date); finput is the input file for raw data
    ID,Month,Year=[],[],[]
    ColDate,ID,Seq= [],[],[]
    ID2,Seq2 = [],[]
    ColDateX,IDX,SeqX= [],[],[]
    ID2X,Seq2X = [],[]
    
    Logfile = []
    
    report1 = "Intial number from Raw Data "+str(len(fnew)/2)+" Seq"
    print(report1) # how many input we have
    b = 0

    for i in np.arange(0,(len(fnew)),2):
        S=fnew[i+1]
        if ('X' not in S) and (len(S)>=1256): #Taken Length above 1256, only bottom limit
            ID.append(fnew[i])
            Seq.append(S)
            w = fnew[i].split("|")
            ColDate.append(w[2]) #collection Date
            b +=1 # b = total number of seq that do not contain X mutation and in this range (1255,1274)
        else:
            IDX.append(fnew[i])
            SeqX.append(S)
            w = fnew[i].split("|")
            ColDateX.append(w[2]) #collection Date

    report2 = "Number of Seq containing X amino acid = "+str(len(fnew)/2 - b)+" Seq"
    report3 = "Free X mutation seq = "+str(b)+" Seq"
    print(report2)
    print(report3)
    print('\n')
    c = 0
    cx = 0
    for i in range(len(ColDate)):#Removing non-specific collection date for Seq no X amino acid
        w = ColDate[i].split("-")
        Month = w[1];Year = w[0]
        if Month != "00" and int(Year) >= 2019:
            ID2.append(ID[i])
            Seq2.append(Seq[i])
            c +=1 #c = total number of seq that show specific month's collection
    
    for i in range(len(ColDateX)):#Removing non-specific collection date for Seq with X amino acid
        wx = ColDateX[i].split("-")
        Month = wx[1];Year = wx[0]
        if Month != "00" and int(Year) >= 2019:
            ID2X.append(ID[i])
            Seq2X.append(Seq[i])
            cx +=1 #c = total number of seq that show specific month's collection
    
    report4 = "Specific collection date seq - Free X = "+str(c)+" Seq"
    report5 = "Non-Specific collection date seq - Free X = "+str(b-c)+" Seq"
    print(report4)
    print(report5)
    print('\n')
    
    report6 = "Specific collection date seq - with X = "+str(cx)+" Seq"
    report7 = "Non-Specific collection date seq - with X = "+str(len(fnew)/2 - b - cx)+" Seq"
    print(report6)
    print(report7)
    print('\n')
    
    Clean =[]
    CleanX =[]
    for i in range(len(ID2)):
        Clean.append(ID2[i])
        Clean.append("\n")
        Clean.append(Seq2[i])
        Clean.append("\n")
        
    for i in range(len(ID2X)):
        CleanX.append(ID2X[i])
        CleanX.append("\n")
        CleanX.append(Seq2X[i])
        CleanX.append("\n")
    
    Logfile.append(report1)
    Logfile.append(report2)
    Logfile.append(report3)
    Logfile.append("\n")
    Logfile.append(report4)
    Logfile.append(report5)
    Logfile.append("\n")
    Logfile.append(report6)
    Logfile.append(report7)
    Logfile.append("\n")
    
    with open('Updated-spikeprot0209-NoX.fasta','w') as f1:#rename this file ----------------------------------------------------------
        for item in Clean:
            f1.write("%s" % item)
    
    with open('Updated-spikeprot0209-withX.fasta','w') as f1X:#rename this file ----------------------------------------------------------
        for item in CleanX:
            f1X.write("%s" % item)
    
    with open('ReportFile.fasta','w') as f2:#rename this file ----------------------------------------------------------
        for item in Logfile:
            f2.write("%s" % item)
            
    return
    
#reference.txt, newdata.fasta,excel.xlsx
def updateunique(reference,newdata,excel): #Obtain an Excel File that contain Previous Unique Sequences with New Freq
    ##############################################################################################
    #Main Backbone of updateunique, NewSeq Function    
    file_input = reference #"Reference0Prot.txt"
    g = open(file_input,'r')
    f0 = g.readlines()
    ref0 = strip_n(f0)
    ref0 = ref0[0]
    
    ProtLen = 1273
    name = newdata[len(newdata)-10:len(newdata)-6] #just for the naming
    
    file_input1 = newdata #"Updated-spikeprot0927.fasta" #New File #read clean file
    g = open(file_input1,'r')
    fnew = g.readlines()
    fnew = strip_n(fnew)
    print("length of the input file = ", len(fnew)/2," Seq")
    print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')
    
    ID,Seq,ColDate,Date,Month,Year,EPI_ID= [],[],[],[],[],[],[]
    YearMon = []
    ColDate2,EPI_ID2,Seq2,ID2 = [],[],[],[]
    label1 = []
    
    
    for i in np.arange(0,(len(fnew)),2):#prepare the ID and ColDate for each data, i for the Full ID, i+1 for Seq
        S=fnew[i+1]
        Seq2.append(S)
        
        w = fnew[i].split("|")
        ColDate.append(w[2])
        EPI_ID2.append(w[3])
        
    for i in range(len(ColDate)):#change non-specific month to 28
        v = ColDate[i].split("-")
        Date = v[2];Month = v[1];Year = v[0]
        if Date == "00":
            Date = "01"
        Collec = Year+"-"+Month+"-"+Date   #Dpt Seq, Epi.ID , Collec date + Mont Year
        YearMon.append(Year+"-"+Month)
        ColDate2.append(Collec)
    
    
    d1 = {'ColDate': ColDate2,'YearMon': YearMon,'EPI_ID': EPI_ID2,'Seq': Seq2, \
              };  df1 = pd.DataFrame(d1) #df new data
    
    
    df1["ColDate"] = pd.to_datetime(df1["ColDate"])#Change the date format str to date
    df1.sort_values(by='ColDate', inplace=True)# time sorted data
    
    require_cols = [1,3,5,6]#Columns that you want to use
    #print (df1)
    
    dfRefMSA = pd.read_excel(excel, sheet_name ='Ref2709', header = 0, \
                        usecols = require_cols)#Reference of unique seq from old data
        #Testing.xlsx
    
    ################################################################################################
    dfUpdateData = pd.merge(df1,dfRefMSA,on="Seq",how="left")#checking a new redundant sequences 
    # if there is some difference, then it is added to the excel
    dfUpdateData = dfUpdateData.fillna(0) #Replace all NaN elements with 0s.
    dfUpdateData = dfUpdateData.drop_duplicates(subset='EPI_ID', keep="first") # Take the first data/earliest seq appears
    ##### count the repeated number of copies in the new data
    MutList= dfUpdateData["MutList"].tolist() #.tolist convert to list
    MutList = [str(item) for item in MutList]     #ini bikin apa sih?
    MutList1, MutListFreq= count_dups(MutList)
    dcount= {'MutList': MutList1,'Freq': MutListFreq \
              };  dfcount = pd.DataFrame(dcount)
    
    
    dfUpdateData = pd.merge(dfUpdateData,dfcount,on="MutList",how="left")#checking a new redundant sequences
    dfUpdateData.sort_values(by='ColDate', inplace=True)# time sorted data
    dfUpdateData = dfUpdateData.drop_duplicates(subset='MutList', keep="first")
    dfUpdateData.to_excel('All-Unique_GISAID-'+'Update-'+name+'.xlsx') #obtain all unique and effective seq + new data
    
    
    print("Number of Unique Seq = ",len(dfUpdateData))
    print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')

def NewSeq(reference,newdata,excel): #Obtain a Fasta Fole that contains the NEW Unique Seq 
    ####################################################################################################    
    #Main Backbone of updateunique, NewSeq Function    
    file_input = reference
    g = open(file_input,'r')
    f0 = g.readlines()
    ref0 = strip_n(f0)
    ref0 = ref0[0]

    ProtLen = 1273

    file_input1 = newdata #New File #read clean file
    #name = newdata[len(newdata)-10:len(newdata)-6]
    g = open(file_input1,'r')
    fnew = g.readlines()
    fnew = strip_n(fnew)
    print("length of the input file = ", len(fnew)/2," Seq")
    print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')
    
    name = newdata[len(newdata)-10:len(newdata)-6] #just for the naming

    ID,Seq,ColDate,Date,Month,Year,EPI_ID= [],[],[],[],[],[],[]
    YearMon = []
    ColDate2,EPI_ID2,Seq2,ID2 = [],[],[],[]
    label1 = []

    for i in np.arange(0,(len(fnew)),2):#prepare the ID and ColDate for each data, i for the Full ID, i+1 for Seq
        S=fnew[i+1]
        Seq2.append(S)
        
        w = fnew[i].split("|")
        ColDate.append(w[2])
        EPI_ID2.append(w[3])
        
    for i in range(len(ColDate)):#change non-specific month to 28
        v = ColDate[i].split("-")
        Date = v[2];Month = v[1];Year = v[0]
        if Date == "00":
            Date = "01"
        Collec = Year+"-"+Month+"-"+Date   #Dpt Seq, Epi.ID , Collec date + Mont Year
        YearMon.append(Year+"-"+Month)
        ColDate2.append(Collec)


    d1 = {'ColDate': ColDate2,'YearMon': YearMon,'EPI_ID': EPI_ID2,'Seq': Seq2, \
              };  df1 = pd.DataFrame(d1) #df new data


    df1["ColDate"] = pd.to_datetime(df1["ColDate"])#Change the date format str to date
    df1.sort_values(by='ColDate', inplace=True)# time sorted data

    require_cols = [1,3,5,6]#Columns that you want to use
    #print (df1)

    dfRefMSA = pd.read_excel(excel, sheet_name ='Ref2709', header = 0, \
                        usecols = require_cols)#Reference of unique seq from old data
        
    ######################################################################################################
    dfUpdateData = pd.merge(df1,dfRefMSA,on="Seq",how="left")#checking a new redundant sequences 
    # if there is some difference, then it is added to the excel
    dfUpdateData = dfUpdateData.fillna(0) #Replace all NaN elements with 0s.
    dfUpdateData = dfUpdateData.drop_duplicates(subset='EPI_ID', keep="first") # Take the first data EPI_ID, throw out duplicated EPI_iD
    
    dfNewUnique = dfUpdateData.loc[dfUpdateData["SeqMSA"]== 0] #Obtain New Seqs
    
    CountNewReport = len(dfNewUnique)
    print("New Report Data: ",CountNewReport)

    #############################################################################
    dfNewUnique["ColDate"] = pd.to_datetime(dfNewUnique["ColDate"])#Change the date format str to date
    dfNewUnique.sort_values(by='ColDate', inplace=True)# time sorted data
    
    dfNewUnique = dfNewUnique.drop_duplicates(subset='Seq', keep="first") #Obtain the UNIQUE New Seq
    
    CountNewUniqueReport = len(dfNewUnique)
    print("New Unique Report Data: ",CountNewUniqueReport)  
    
    dfNewUnique.to_excel('Processed-NewUnique-'+name+'.xlsx')
    
    dfEPI = dfNewUnique['EPI_ID']
    dfSeq = dfNewUnique['Seq']
    lEPI = dfEPI.values.tolist()
    lSeq = dfSeq.values.tolist()
    
    UniqueSeqTXT = []
    MSAref = 'MSARef.txt'
    reflist = Split(MSAref)
    
    for i in range(len(reflist)):
        UniqueSeqTXT.append(reflist[i])
    
    for i in range(len(lEPI)):
        UniqueSeqTXT.append('>'+lEPI[i])
        UniqueSeqTXT.append('\n')
        UniqueSeqTXT.append(lSeq[i])
        UniqueSeqTXT.append('\n')
    
    with open('ProcessedNewUnique-'+name+'.fasta','w') as f1:#rename this file ----------------------------------------------------------
        for item in UniqueSeqTXT:
            f1.write("%s" % item)
    ########################################################################################
    print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')

def Read_AllignmentMSA(fastafile): #Obtain a txt file that has ID and Seq of the MSA
    file_input = fastafile#Input File after MSA

    g = open(file_input,'r')

    f1 = g.readlines()
    length = len(f1) #number of rows in fasta file

    cnt,add, cnte= [],[],[]
    add.append(f1[0])
    b = 0 #initial row which contain ">"
    for i in range(1,length) :
        for l in f1[i]:
            if l[0] == '>': #if f1[i] contain ">"
                add.append(f1[i])
                f1[i] =""            
                q = i
                cnt = f1[b+1:q]
                for element in cnt:
                    cnte.append(element.strip())
                cnte.append("\n")
                b = q
    cnt = f1 [b+1:length]
    for element in cnt:
        cnte.append(element.strip())
    #cnte.append("\n")

    with open('Read_Alignment-Seq.txt', 'w') as f:
    	for item in cnte:
    		f.write("%s" % item)

    with open('Read_Alignment-ID.txt', 'w') as f:
    	for item in add:
    		f.write("%s" % item)

    ID_input = 'Read_Alignment-ID.txt' #Contain only ID
    h = open(ID_input,'r')
    fID = h.readlines()
    fnew = strip_n(fID)

    
    newfile_input = 'Read_Alignment-Seq.txt' # Contain only Seq
    g = open(newfile_input,'r')
    fnew = g.readlines()
    fnew = strip_n(fnew)
    print(len(fID))
    print(len(fnew))
    
    ID_and_Seq = []
    for i in range(len(fID)):
        ID_and_Seq.append(fID[i])
        ID_and_Seq.append("\n")
        ID_and_Seq.append(fnew[i])
        ID_and_Seq.append("\n")

    """ THE FINAL OUTPUT """ #Contain both Seq and ID txt
    with open('No_Up_lim_After_MSA-Seq_All.txt', 'w') as f:#change the number based on the txt file
        for item in ID_and_Seq:
            f.write("%s" % item)

def Read_AllignmentExcel(file): #Convert the Read_AllingmentMSA txt file to excel
    
    file_input = file #Input File after MSA
    g = open(file_input,'r')
    f1 = g.readlines()
    
    LenSeq = []
    ID,Seq = ReadAlign(f1)
    
    for i in range(len(Seq)):
        LenSeq.append(len(Seq[i]))
    
    d1 = {'ID': ID,'Seq': Seq,'Length': LenSeq};df1 = pd.DataFrame(d1)
    d2 = {'Seq': Seq};df2 = pd.DataFrame(d2)
    df2 = SplitChar(df2)
    df = pd.concat([df1, df2], axis=1)
    
    df.to_excel('NoUpLimMSAJanSeq_Trial3.xlsx')#rename it
    
def UpdatingNewData(PreProcessingResult, testingexcel):

    file_input1 = PreProcessingResult #"Updated-spikeprot0927.fasta" #New File #read clean file
    g = open(file_input1,'r')
    fnew = g.readlines()
    fnew = strip_n(fnew)
    print("length of the input file = ", len(fnew)/2," Seq")
    print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')
    
    name = file_input1[len(file_input1)-10:len(file_input1)-6] #just for the naming

    ID,Seq,ColDate,Date,Month,Year,EPI_ID= [],[],[],[],[],[],[]
    YearMon = []
    ColDate2,EPI_ID2,Seq2,ID2 = [],[],[],[]
    label1 = []
    
    for i in np.arange(0,(len(fnew)),2):#prepare the ID and ColDate for each data, i for the Full ID, i+1 for Seq
        S=fnew[i+1]
        Seq2.append(S)
        
        w = fnew[i].split("|")
        ColDate.append(w[2])
        EPI_ID2.append(w[3])
        
    for i in range(len(ColDate)):#change non-specific month to 28
        v = ColDate[i].split("-")
        Date = v[2];Month = v[1];Year = v[0]
        if Date == "00":
            Date = "01"
        Collec = Year+"-"+Month+"-"+Date   #Dpt Seq, Epi.ID , Collec date + Mont Year
        YearMon.append(Year+"-"+Month)
        ColDate2.append(Collec)
    
    
    d1 = {'ColDate': ColDate2,'YearMon': YearMon,'EPI_ID': EPI_ID2,'Seq': Seq2, \
              };  df1 = pd.DataFrame(d1) #df new data
    
    
    df1["ColDate"] = pd.to_datetime(df1["ColDate"])#Change the date format str to date
    df1.sort_values(by='ColDate', inplace=True)# time sorted data
    
    
    require_cols = [1,3,5,6]#Columns that you want to use
    #print (df1)
    
    dfRefMSA = pd.read_excel(testingexcel, sheet_name ='Ref2709', header = 0, \
                        usecols = require_cols)#Reference of unique seq from old data
        
    dfUpdateData = pd.merge(df1,dfRefMSA,on="Seq",how="left")
    # if there is some difference, then it is added to the excel
    dfUpdateData = dfUpdateData.fillna(0) #Replace all NaN elements with 0s.
    dfUpdateData = dfUpdateData.drop_duplicates(subset='EPI_ID', keep="first")
    listSeq0 = dfUpdateData['Seq'].values.tolist()
    listSeq1, listfreque = count_dups(listSeq0)
    dcount= {'Seq': listSeq1,'Freq': listfreque \
              };  dfcount = pd.DataFrame(dcount)
        
    dfUpdateData = pd.merge(dfUpdateData,dfcount,on="Seq",how="left")
    dfUpdateData = dfUpdateData.fillna(0) #Replace all NaN elements with 0s.
    dfUpdateData = dfUpdateData.drop_duplicates(subset='EPI_ID', keep="first")
    dfUpdateData.sort_values(by='ColDate', inplace=True)
    dfUpdateData = dfUpdateData.drop_duplicates(subset='Seq', keep="first")
    print(len(dfUpdateData))
    s = 0
    length = int(len(dfRefMSA))
    gap = int(len(dfRefMSA)/3)


    for i in range(0,length+1,gap):
        print(str(i), "-Iteration")
        dfprint = dfUpdateData[i:(i+gap)]
        dfprint.to_excel('All_Unique_Freque-RETRY'+str(i)+'-'+name+'.xlsx') #obtain all unique and effective seq + new data
        
def MonthlyData(PreProcessingResult, testingexcel):
    file_input1 = PreProcessingResult #"Updated-spikeprot0927.fasta" #New File #read clean file
    g = open(file_input1,'r')
    fnew = g.readlines()
    fnew = strip_n(fnew)
    print("length of the input file = ", len(fnew)/2," Seq")
    print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')
    
    
    name = file_input1[len(file_input1)-10:len(file_input1)-6] #just for the naming

    ID,Seq,ColDate,Date,Month,Year,EPI_ID= [],[],[],[],[],[],[]
    YearMon = []
    ColDate2,EPI_ID2,Seq2,ID2 = [],[],[],[]
    label1 = []
    
    for i in np.arange(0,(len(fnew)),2):#prepare the ID and ColDate for each data, i for the Full ID, i+1 for Seq
        S=fnew[i+1]
        Seq2.append(S)
        
        w = fnew[i].split("|")
        ColDate.append(w[2])
        EPI_ID2.append(w[3])
        
    for i in range(len(ColDate)):#change non-specific month to 28
        v = ColDate[i].split("-")
        Date = v[2];Month = v[1];Year = v[0]
        if Date == "00":
            Date = "01"
        Collec = Year+"-"+Month+"-"+Date   #Dpt Seq, Epi.ID , Collec date + Mont Year
        YearMon.append(Year+"-"+Month)
        ColDate2.append(Collec)
    
    
    d1 = {'ColDate': ColDate2,'YearMon': YearMon,'EPI_ID': EPI_ID2,'Seq': Seq2, \
              };  df1 = pd.DataFrame(d1) #df new data
    
    
    df1["ColDate"] = pd.to_datetime(df1["ColDate"])#Change the date format str to date
    df1.sort_values(by='ColDate', inplace=True)# time sorted data
    
    require_cols = [1,3,5,6]#Columns that you want to use
    #print (df1)
    
    dfRefMSA = pd.read_excel(testingexcel, sheet_name ='Ref2709', header = 0, \
                        usecols = require_cols)#Reference of unique seq from old data
        
    dfUpdateData = pd.merge(df1,dfRefMSA,on="Seq",how="left")
    # if there is some difference, then it is added to the excel
    dfUpdateData = dfUpdateData.fillna(0) #Replace all NaN elements with 0s.
    dfUpdateData = dfUpdateData.drop_duplicates(subset='EPI_ID', keep="first")
    YearMon1, YearMonFreq = count_dups(YearMon)#obtain the effective month(Dec'19 - Sept'21) 
    d= {'YearMon': YearMon1,'Freq': YearMonFreq \
              };  dfYearMon = pd.DataFrame(d)
    # print(dfYearMon)
    
    dfPassHor = pd.DataFrame(columns = ['ColDate', 'YearMon', 'EPI_ID', 'Seq', 'SeqMSA', 'N.mut', 'MutList', 'Freq'])
    print(dfPassHor)
    lpass =  0
    for i in range(len(YearMon1)):
    
        select = YearMon1[i]
    
        dfMonthly = dfUpdateData.loc[dfUpdateData["YearMon"] == select]
        
        MutList = dfMonthly["MutList"].tolist()
        MutList = [str(item) for item in MutList]
        MutList1, MutListFreq= count_dups(MutList)
        dcount= {'MutList': MutList1,'Freq': MutListFreq \
              };  dfcount = pd.DataFrame(dcount)
      
        
        dfprint = pd.merge(dfMonthly,dfcount,on="MutList",how="left")
        dfprint["ColDate"] = pd.to_datetime(dfprint["ColDate"])#Change the date format str to date
        dfprint.sort_values(by='ColDate', inplace=True)# time sorted data
        
        dfprint = dfprint.drop_duplicates(subset='MutList', keep="first")#Keep the eff seq
        lpass += len(dfprint)
        print("iteration -",select, "Number of Seq = ", len(dfprint))# "|Cumulative =",lpass
    
        # dfFinal = pd.concat([dfprint, df3, df4], axis=1)# concatinating 3 df
        # dfprint.to_excel('Retrieve_GISAID-Update0927-'+select+'.xlsx')
        dfPassHor["space"] = ""
        dfPassHor = pd.concat([dfPassHor,dfprint], axis=0, join="outer")
    print(len(dfPassHor))
    #dfPassHor.to_excel('All_Monthly-Unique-Update+'+name+'.xlsx')
    
    print('\n',lpass)
    s = 0
    length = int(len(dfPassHor))
    gap = int(len(dfPassHor)/3)

    for i in range(0,length+1,gap):
        print(str(i), "-Iteration")
        dfprint = dfPassHor[i:(i+gap)]
        dfprint.to_excel('Monthly_Unique_Freque'+str(i)+'-'+name+'.xlsx')        

def UniqueSeq(reference,newdata,name): #Obtain an Excel File that contain Previous Unique Sequences with New Freq
    ##############################################################################################
    #Main Backbone of updateunique, NewSeq Function    
    ref0 = strip_n( open(reference,'r').readlines() ) #"Reference0Prot.txt"
    ref0 = ref0[0]
    
    # name = newdata[len(newdata)-10:len(newdata)-6] #just for the naming
     
    fnew = strip_n( open(newdata,'r').readlines() ) #"Updated-spikeprot0927.fasta" #New File #read clean file
    print("length of the input file = ", len(fnew)/2," Seq")
    print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')
    
    ID,Seq,ColDate,Date,Month,Year,EPI_ID= [],[],[],[],[],[],[]
    YearMon = []
    ColDate2,EPI_ID2,Seq2,ID2 = [],[],[],[]
    
    
    for i in np.arange(0,(len(fnew)),2):#prepare the ID and ColDate for each data, i for the Full ID, i+1 for Seq
        S=fnew[i+1]
        Seq2.append(S)
        
        w = fnew[i].split("|")
        ColDate.append(w[2])
        EPI_ID2.append(w[3])
        
    for i in range(len(ColDate)):#change non-specific month to 28
        v = ColDate[i].split("-")
        Date = v[2];Month = v[1];Year = v[0]
        if Date == "00":
            Date = "01"
        Collec = Year+"-"+Month+"-"+Date   #Dpt Seq, Epi.ID , Collec date + Mont Year
        YearMon.append(Year+"-"+Month)
        ColDate2.append(Collec)
    
    
    d1 = {'ColDate': ColDate2,'YearMon': YearMon,'EPI_ID': EPI_ID2,'Seq': Seq2, \
              };  df1 = pd.DataFrame(d1) #df new data
    
    
    df1["ColDate"] = pd.to_datetime(df1["ColDate"])#Change the date format str to date
    df1.sort_values(by='ColDate', inplace=True)# time sorted data
    
    ################################################################################################

    dfUpdateData = df1.drop_duplicates(subset='Seq', keep="first") # Take the first data/earliest seq appears
    # ##### count the repeated number of copies in the new data
    # Seq= dfUpdateData["Seq"].tolist() #.tolist convert to list
    # # Seq = [str(item) for item in Seq]     #ini bikin apa sih?
    # Seq1, SeqFreq= count_dups(Seq)
    # dcount= {'Seq': Seq1,'Freq': Seq};  dfcount = pd.DataFrame(dcount)
    
    # dfUpdateData = pd.merge(dfUpdateData,dfcount,on="Seq",how="left")#checking a new redundant sequences
    # dfUpdateData.sort_values(by='ColDate', inplace=True)# time sorted data
    # dfUpdateData = dfUpdateData.drop_duplicates(subset='Seq', keep="first")
    dfUpdateData.to_excel('All-Unique_GISAID-'+'Update-'+name+'.xlsx') #obtain all unique and effective seq + new data
    
    print("Number of Unique Seq = ",len(dfUpdateData))
    print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')
    
print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')
ProtLen = 1273   
####################################################################
# December 1218 Data
# MonthlyData('Updated-spikeprot1218.fasta','UpdatedTesting.xlsx')
# file_input = 'spikeprot1218.fasta'
# g = open(file_input,'r')
# f0 = g.readlines()
# f0 = strip_n(f0)
# PreProcessing(f0)


# updateunique('Reference0Prot.txt','Updated-spikeprot0927.fasta','Testing.xlsx')
# NewSeq('Reference0Prot.txt','Updated-spikeprot1218.fasta','Testing.xlsx')
# Separate('ProcessedNewUnique-1218.fasta',100)
# ReadOutFile('slurm-548417.out')
# Read_AllignmentMSA('msa1218-x.fa')
# SearchDel('After_MSA-Seq1218_100.txt')
# SearchMoreThan2Del('After_MSA-Seq1218_100.txt')
# Read_AllignmentExcel('testingdata_o.fa')
# Read_AllignmentMSA('msa1218_allSeq.fa')

# Separate('Modified_CompiledAfterMSA.txt',10000)

# UpdatingNewData('Updated-spikeprot1218.fasta','CheckEverythingUpdateTesting.xlsx')
# MonthlyData('Updated-spikeprot1218.fasta','CheckEverythingUpdateTesting.xlsx')

####################################################################################
# February 0209 Data
f0 = strip_n( open('spikeprot0209.fasta','r').readlines() )
# PreProcessing(f0)
UniqueSeq('Reference0Prot.txt','Updated-spikeprot0209-NoX.fasta','0209')
# NewSeq('Reference0Prot.txt','No_Up_Lim_Updated-spikeprot0209.fasta','CheckEverythingUpdateTesting.xlsx')

########### Place NewSeq input to MSA ##################################################
# Read_AllignmentMSA('no_limit_msa0103_allSeq.fa')
####At this step, if you fail the code, you can split the .fa file to many files, use Separate() Function

# Read_AllignmentExcel('No_Up_lim_After_MSA-Seq_All.txt')
######################################################################################################

# UpdatingNewData('No_Up_Lim_Updated-spikeprot0103.fasta','0103UpdateTesting.xlsx')
# MonthlyData('No_Up_Lim_Updated-spikeprot0103.fasta','0103UpdateTesting.xlsx')

########################################################################
# December 0204 
# file_input = 'spikeprot0204.fasta'
# g = open(file_input,'r')
# f0 = g.readlines()
# f0 = strip_n(f0)
# PreProcessing(f0)

# updateunique('Reference0Prot.txt','Updated-spikeprot0204.fasta','0103UpdateTesting.xlsx')
# NewSeq('Reference0Prot.txt','Updated-spikeprot0204.fasta','0103UpdateTesting.xlsx')


print('\nTIME TAKEN: ' + str(time.process_time() - start_time) + 's')
