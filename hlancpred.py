###################################################################################
# HLAncPred is developed for predicting and scanning class-I non-classical HLA      #
# binder peptides. It is developed by Prof G. P. S. Raghava's group.                #
# Please cite: HLAncPred; available at https://webs.iiitd.edu.in/raghava/hlancpred/ #
##################################################################################
import argparse
import warnings
import subprocess
import pkg_resources
import os
import sys
import numpy as np
import pandas as pd
import math
import itertools
from collections import Counter
import pickle
import re
import glob
import time
from time import sleep
from tqdm import tqdm
import xgboost
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='Please provide following arguments') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code")
parser.add_argument("-a", "--allele", type=str.upper, choices = ['G0101','G0103','G0104','E0101','E0103'],required=True, help="Please provide the name of allele for the prediction of binder peptides")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
parser.add_argument("-j", "--job",type=int, choices = [1,2], help="Job Type: 1:predict, and 2:scan, by default 1")
parser.add_argument("-w","--winleng", type=int, choices =range(8, 16), help="Window Length: 8 to 35 (scan mode only), by default 9")
parser.add_argument("-d","--display", type=int, choices = [1,2], help="Display: 1:Only binders, 2: All peptides, by default 1")
args = parser.parse_args()

# Function for generating pattern of a given length
def seq_pattern(aa_seq,win_len):
    aa_seq == aa_seq.upper
    seq_pat=[]
    for i1 in range(0, (len(aa_seq) + 1 - win_len)):
        i2 = i1 + int(win_len)
        seq_pat += [aa_seq[i1:i2]]
    return seq_pat

#Feature Generation
def feature_gen(file):
    def equal_length(file):
        df2 = pd.read_csv(file,header=None,names=['Seq'])
        cc = []
        for i in range(0,len(df2)):
            cc.append(df2['Seq'][i]+(15-len(df2['Seq'][i]))*'X')
        df3 = pd.DataFrame(cc)
        return df3
    def aab(file):
        std = list("ACDEFGHIKLMNPQRSTVWYX")
        filename, file_extension = os.path.splitext(file)
        df = equal_length(file)
        uu = []
        for ss in df[0]:
            uu.append(len(ss))
        zz = df.iloc[:,0]
        f = open('sam_allcomp.aab', mode='w')
        sys.stdout = f
        A=('1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
        C=('0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
        D=('0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
        E=('0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
        F=('0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
        G=('0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
        H=('0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
        I=('0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0')
        K=('0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0')
        L=('0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0')
        M=('0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0')
        N=('0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0')
        P=('0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0')
        Q=('0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0')
        R=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0')
        S=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0')
        T=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0')
        V=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0')
        W=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0')
        Y=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0')
        X=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1')
        for mm in range (1,max(uu)+1):
            for ee in std:
                print(ee+str(mm),end=',')
        print("")
        for i in range(0,len(zz)):
            for j in zz[i]:
                if j == "A":
                    print(''.join(A), end = ',')
                if j == "C":
                    print(''.join(C), end = ',')
                if j == "D":
                    print(''.join(D), end = ',')
                if j == "E":
                    print(''.join(E), end = ',')
                if j == "F":
                    print(''.join(F), end = ',')
                if j == "G":
                    print(''.join(G), end = ',')
                if j == "H":
                    print(''.join(H), end = ',')
                if j == "I":
                    print(''.join(I), end = ',')
                if j == "K":
                    print(''.join(K), end = ',')
                if j == "L":
                    print(''.join(L), end = ',')
                if j == "M":
                    print(''.join(M), end = ',')
                if j == "N":
                    print(''.join(N), end = ',')
                if j == "P":
                    print(''.join(P), end = ',')
                if j == "Q":
                    print(''.join(Q), end = ',')
                if j == "R":
                    print(''.join(R), end = ',')
                if j == "S":
                    print(''.join(S), end = ',')
                if j == "T":
                    print(''.join(T), end = ',')
                if j == "V":
                    print(''.join(V), end = ',')
                if j == "W":
                    print(''.join(W), end = ',')
                if j == "Y":
                    print(''.join(Y), end = ',')
                if j == "X":
                    print(''.join(X), end = ',')
            print("")
        f.truncate()
    aab(file)
    df1 = pd.read_csv("sam_allcomp.aab")
    df19 = df1.iloc[:,:-1]
    filelist=glob.glob("sam_allcomp*")
    for file_2 in filelist:
        os.remove(file_2)
    return df19

def adjusted_classes(y_scores, t):
    return [1 if y >= t else 0 for y in y_scores]

def Perform_testing(clf,name,X,t):
    Y_pred = clf.predict(X)
    Y_scores=[]
    Y_scores=clf.predict_proba(X)[:,-1]
    Y_pred = adjusted_classes(Y_scores,t)
    return Y_pred,Y_scores

print('####################################################################################')
print('# This program HLAncPred is developed for predicting and scanning non-classical    #')
print('# class-I HLA binder peptides, developed by Prof G. P. S. Raghava group.           #')
print('# Please cite: HLAncPred; available at https://webs.iiitd.edu.in/raghava/hlancpred/  #')
print('####################################################################################')

# Parameter initialization or assigning variable for command level arguments

Sequence= args.input        # Input variable 
al = args.allele
# Output file 
 
if args.output == None:
    result_filename= "outfile.csv" 
else:
    result_filename = args.output
         
# Job Type 
if args.job == None:
        Job = int(1)
else:
        Job = int(args.job)
# Window Length 
if args.winleng == None:
        Win_len = int(9)
else:
        Win_len = int(args.winleng)

# Display
if args.display == None:
        dplay = int(1)
else:
        dplay = int(args.display)


if Job==2:
    print("\n");
    print('##############################################################################')
    print('Summary of Parameters:')
    print('Input File: ',Sequence,'Chosen Allele: ', al,'; Job Type: ',Job)
    print('Output File: ',result_filename,'; Window Length: ',Win_len,'; Display: ',dplay)
    print('##############################################################################')
else:
    print("\n");
    print('##############################################################################')
    print('Summary of Parameters:')
    print('Input File: ',Sequence,'; Chosen Allele: ', al,'; Job Type: ',Job)
    print('Output File: ',result_filename,'; Display: ',dplay)
    print('# ############################################################################')
#------------------ Read input file ---------------------
def load_model(path):
    clf = pickle.load(open(path,'rb'))
    return clf

f=open(Sequence,"r")
len1 = f.read().count('>')
f.close()

f=open(Sequence,"r")
seqs=[]
seqid=[]
str1='';
header=0
line=0
if len1 >= 1: # read fasta file
    for l in f:
        if l.startswith('>'):
            if header != 0:
                seqs += [str1]
                str1 = '';
            header = 1
            line+=1
            seqid += [l.rstrip()]
        else:
            str1 += l.rstrip()
    seqs += [str1]
else: # read single line file    
    for l in f:
        if len(l) >= 8 and len(l)<= 15:
            seqs+=[l.strip()]
            seqid += ['>Seq_' + str(line)]
        line+=1
f.close()

fout= open(result_filename,"w+")

i1 = 0
#======================= Prediction Module start from here =====================
if Job == 1:
    print('\n======= Thanks for using Predict module of HLAncPred. Your results will be stored in file :',result_filename,' =====\n')
    fout.write('# Sequence_ID, Sequence, Prediction, Score\n')
    for Sequence in tqdm(seqs):
        header = seqid[i1]
        if len(Sequence) >= 8: 
            if len(Sequence) >= 15:
                Sequence = Sequence[0:15]
            ss = []
            ss.append(Sequence)
            pd.DataFrame(ss).to_csv(Sequence,index=None,header=False)
            X = feature_gen(Sequence)
            if al == 'G0101':
                clf=load_model('Models/SVC_G0101.pkl')
                Y_pred,Y_score=Perform_testing(clf,'SVC',X,0.43)
                Y_score = np.round(Y_score,2)
                flag=""
                if Y_pred[0]==1:
                    flag='HLA-G*01:01 Binder'
                else:
                    flag='HLA-G*01:01 Non-binder'
                if dplay == 1:
                    if Y_pred[0]==1:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                else:
                    fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
            if al == 'G0103':
                clf=load_model('Models/XGB_G0103.pkl')
                Y_pred,Y_score=Perform_testing(clf,'XGB',X,0.53)
                Y_score = np.round(Y_score,2)
                flag=""
                if Y_pred[0]==1:
                    flag='HLA-G*01:03 Binder'
                else:
                    flag='HLA-G*01:03 Non-binder'
                if dplay == 1:
                    if Y_pred[0]==1:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                else:
                    fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
            if al == 'G0104':
                clf=load_model('Models/ET_G0104.pkl')
                Y_pred,Y_score=Perform_testing(clf,'ET',X,0.44)
                Y_score = np.round(Y_score,2)
                flag=""
                if Y_pred[0]==1:
                    flag='HLA-G*01:04 Binder'
                else:
                    flag='HLA-G*01:04 Non-binder'
                if dplay == 1:
                    if Y_pred[0]==1:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                else:
                    fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
            if al == 'E0101':
                clf=load_model('Models/ET_E0101.pkl')
                Y_pred,Y_score=Perform_testing(clf,'ET',X,0.5)
                Y_score = np.round(Y_score,2)
                flag=""
                if Y_pred[0]==1:
                    flag='HLA-E*01:01 Binder'
                else:
                    flag='HLA-E*01:01 Non-binder'
                if dplay == 1:
                    if Y_pred[0]==1:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                else:
                    fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
            if al == 'E0103':
                clf=load_model('Models/RF_E0103.pkl')
                Y_pred,Y_score=Perform_testing(clf,'RF',X,0.49)
                Y_score = np.round(Y_score,2)
                flag=""
                if Y_pred[0]==1:
                    flag='HLA-E*01:03 Binder'
                else:
                    flag='HLA-E*01:03 Non-binder'
                if dplay == 1:
                    if Y_pred[0]==1:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                else:
                    fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
            os.remove(Sequence)
            i1 = i1 +1
#=============== Scan Model start from here ==================
else:
    print('\n======= Thanks for using Scan module of HLAncPred. Your results will be stored in file :',result_filename,' =====\n')
    print('==== Scanning Peptides: Processing sequences please wait ...')
    i1 = 0
    fout.write('# Sequence_ID, pattern, Prediction, Score\n')
    for Sequence1 in tqdm(seqs):
        pat_seq=[]
        header = seqid[i1]
        if len(Sequence1) >= Win_len: 
            pat_seq = seq_pattern(Sequence1,Win_len)
            for Sequence in pat_seq:
                ss = []
                ss.append(Sequence)
                pd.DataFrame(ss).to_csv(Sequence,index=None,header=False)
                X = feature_gen(Sequence)
                if al == 'G0101':
                    clf=load_model('Models/SVC_G0101.pkl')
                    Y_pred,Y_score=Perform_testing(clf,'SVC',X,0.43)
                    Y_score = np.round(Y_score,2)
                    flag=""
                    if Y_pred[0]==1:
                        flag='HLA-G*01:01 Binder'
                    else:
                        flag='HLA-G*01:01 Non-binder'
                    if dplay == 1:
                        if Y_pred[0]==1:
                            fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                    else:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                if al == 'G0103':
                    clf=load_model('Models/XGB_G0103.pkl')
                    Y_pred,Y_score=Perform_testing(clf,'XGB',X,0.53)
                    Y_score = np.round(Y_score,2)
                    flag=""
                    if Y_pred[0]==1:
                        flag='HLA-G*01:03 Binder'
                    else:
                        flag='HLA-G*01:03 Non-binder'
                    if dplay == 1:
                        if Y_pred[0]==1:
                            fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                    else:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                if al == 'G0104':
                    clf=load_model('Models/ET_G0104.pkl')
                    Y_pred,Y_score=Perform_testing(clf,'ET',X,0.44)
                    Y_score = np.round(Y_score,2)
                    flag=""
                    if Y_pred[0]==1:
                        flag='HLA-G*01:04 Binder'
                    else:
                        flag='HLA-G*01:04 Non-binder'
                    if dplay == 1:
                        if Y_pred[0]==1:
                            fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                    else:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                if al == 'E0101':
                    clf=load_model('Models/ET_E0101.pkl')
                    Y_pred,Y_score=Perform_testing(clf,'ET',X,0.5)
                    Y_score = np.round(Y_score,2)
                    flag=""
                    if Y_pred[0]==1:
                        flag='HLA-E*01:01 Binder'
                    else:
                        flag='HLA-E*01:01 Non-binder'
                    if dplay == 1:
                        if Y_pred[0]==1:
                            fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                    else:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                if al == 'E0103':
                    clf=load_model('Models/RF_E0103.pkl')
                    Y_pred,Y_score=Perform_testing(clf,'RF',X,0.49)
                    Y_score = np.round(Y_score,2)
                    flag=""
                    if Y_pred[0]==1:
                        flag='HLA-E*01:03 Binder'
                    else:
                        flag='HLA-E*01:03 Non-binder'
                    if dplay == 1:
                        if Y_pred[0]==1:
                            fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                    else:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                os.remove(Sequence)
        i1 = i1 +1
    fout.close()
