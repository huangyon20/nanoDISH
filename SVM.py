import os
import argparse
from sklearn.svm import OneClassSVM
import pandas as pd
import numpy as np
from collections import Counter
pd.options.mode.chained_assignment = None


def load_fasta(seqFn):
    """
    seqFn               -- Fasta file
    Return:
        {tid1: seq1, ...} 
    """
    fasta = {}
    cur_tid = ''
    cur_seq = ''
    
    for line in open(seqFn):
        if line[0] == '>':
            if cur_seq != '':
                fasta[cur_tid] = cur_seq
                cur_seq = ''
            data = line[1:].split(None, 1)
            cur_tid = data[0]
        else:
            cur_seq += line.rstrip()
    
    if cur_seq != '':
        fasta[cur_tid] = cur_seq
    return fasta

def SVM(testevent,trainevent,poslist,KN,GA,NU,cols):
    Staticfile={}
    TestwithPred=pd.DataFrame()
    for b in poslist: #对每个坐标，提取出所有reads在该位点的值，判断两个群之间的差异
        Train=trainevent[trainevent['position']==b]
        Test=testevent[testevent['position']==b]
        if Test.shape[0] >10 and Train.shape[0] >10:
            #training of unmodified
            clf = OneClassSVM(kernel=KN,gamma=GA,nu=NU).fit(np.array(Train.iloc[:,cols]))
            #Predicting of modified
            Testout=clf.predict(np.array(Test.iloc[:,cols]))
            Test['Predict']=Testout
            TestwithPred=pd.concat([TestwithPred,Test],axis=0)
            #统计-1和1的个数
            d = Counter(Testout)
            d_s = sorted(d.items(),key=lambda x:x[1],reverse=True)
            #将结果输出到Staticfile字典
            if len(d_s)==2:
                Staticfile[b]=[d_s[1][1],Train['reference_kmer'].tolist()[1][0:1],Test.shape[0],Train.shape[0]]
            if len(d_s)==1:
                Staticfile[b]=[0,Train['reference_kmer'].tolist()[1][0:1],Test.shape[0],Train.shape[0]]
        else:
            Staticfile[b]=[0,'NULL',Test.shape[0],Train.shape[0]]
    return Staticfile,TestwithPred


def args():
    """
    predict modification of each bases by comparing the event file from Mod and Unmod group

    """ 
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--Mod_Event",   type=str, required=True,      help='The Modified collapsed event file from dataprocessing.py and picksize.py.')
    parser.add_argument("-u", "--Unmod_Event", type=str, required=True,      help='The Unmodified collapsed event file from dataprocessing.py and picksize.py.')
    parser.add_argument("-r", "--reference",   type=str, required=True,      help='The reference file.')
    parser.add_argument("-s", "--seq_name",    type=str, required=True,      help='The specific sequence name in the reference file.')
    parser.add_argument("-c", "--columns",     type=list, default=[4,5,6],   help='the features used for SVM, 4: current, 5: current stdv, 6: dwell time. The default is [4,5,6].')
    parser.add_argument("-KN", "--kernel",     type=str, default='rbf',      help='the SVM parameter: kernel ，the default is (rbf).')
    parser.add_argument("-GA", "--gamma",      type=str, default='scale',    help='the SVM parameter: gamma ，the default is (scale).')
    parser.add_argument("-NU", "--nu",         type=float, default=0.01,     help='the SVM parameter: gamma ，the default is (0.01).')
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = args()

    if not os.path.exists(args.Mod_Event):
        print(f'Modified collapsed event file not found: {args.Mod_Event}')
        exit()
    if not os.path.exists(args.Unmod_Event):
        print(f'Unmodified collapsed event file not found: {args.Unmod_Event}')
        exit()
    if not os.path.exists(args.reference):
        print(f'Reference file not found: {args.reference}')
        exit()
    assert args.columns in ([4,5], [4,5,6],[4,6]), "method should be one of [4,5], [4,5,6], [4,6]."

    #get the positions of specific sequence 
    seqs=load_fasta(args.reference)
    length=len(seqs[args.seq_name])
    poslist=[i for i in range(1,length+1)] 
    # read in the Event files
    Mod_Event= pd.read_csv(args.Mod_Event,sep='\t')
    Unmod_Event= pd.read_csv(args.Unmod_Event,sep='\t')
    # predict every bases, add the results to the Mod_event file as a new column
    Staticfile,TestwithPred=SVM(Mod_Event,Unmod_Event,poslist,args.kernel,args.gamma,args.nu,args.columns)
    # write the results to files
    Pred_event=args.Mod_Event+'.Mod_Profile'
    Pred_static=args.Mod_Event+'.static'
    TestwithPred.to_csv(Pred_event,sep='\t',header=True, index=True)
    out=pd.DataFrame.from_dict(Staticfile, orient='index',columns=['Modifired_NO.','Kmer','Test_NO.','Train_NO.'])
    out.to_csv(Pred_static,sep='\t',index=True,)

