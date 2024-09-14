import os
import argparse
from sklearn.svm import OneClassSVM
import pandas as pd
import numpy as np
from collections import Counter
pd.options.mode.chained_assignment = None


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
    parser.add_argument("-m", "--Mod_Event",   type=str, required=True,      help='The Modified collapsed event file.')
    parser.add_argument("-u", "--Unmod_Event", type=str, required=True,      help='The Unmodified collapsed event file.')
    parser.add_argument("-l", "--length",      type=int, required=True,      help='The RNA intermidiate length.')
    parser.add_argument("-b", "--bias",        type=int, default=0,          help='The read length bias. Set the same value with pickintermediate.py, if pickintermediate.py be used.')
    parser.add_argument("-c", "--columns",     type=int, nargs='+', default=[4,5,6],   help='the features used for SVM, 4: current, 5: current stdv, 6: dwell time. The default is [4,5,6].')
    parser.add_argument("-KN", "--kernel",     type=str, default='rbf',      help='The SVM parameter: kernel ，the default is (rbf).')
    parser.add_argument("-GA", "--gamma",      type=float, default=0.01,     help='The SVM parameter: gamma ，the default is (0.01).')
    parser.add_argument("-NU", "--nu",         type=float, default=0.01,     help='The SVM parameter: nu ，the default is (0.01).')
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
    assert args.columns in ([4,5], [4,5,6],[4,6]), "method should be one of [4,5], [4,5,6], [4,6]."

    #get the positions of specific sequence 
    poslist=[i for i in range(1,args.length+args.bias+1)]
    # read in the Event files
    Mod_Event= pd.read_csv(args.Mod_Event,sep='\t',index_col=0)
    Unmod_Event= pd.read_csv(args.Unmod_Event,sep='\t',index_col=0)
    # predict every bases, add the results to the Mod_event file as a new column
    print(args.columns)
    Staticfile,TestwithPred=SVM(Mod_Event,Unmod_Event,poslist,args.kernel,args.gamma,args.nu,args.columns)
    # write the results to files
    Pred_event=args.Mod_Event+'.Mod_Profile'
    Pred_static=args.Mod_Event+'.static'
    TestwithPred.to_csv(Pred_event,sep='\t',header=True, index=True)
    out=pd.DataFrame.from_dict(Staticfile, orient='index',columns=['Modifired_NO.','Kmer','Test_NO.','Train_NO.'])
    out.to_csv(Pred_static,sep='\t',index=True)

