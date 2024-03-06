import os
import argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
pd.options.mode.chained_assignment = None
import scipy.stats
from statistics import mean
from scipy import stats

from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.cluster import DBSCAN

def load_fasta(seqFn):
    """
    seqFn               -- Fasta file
    Return dict:
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

def transPredictedEventToMatrix(TestwithPred,poslist): #SVM步骤中产生的TestwithPred文件含有predict的信号值，1和-1，可以将这些值作为表达矩阵，用于reads按照结构的聚类
    """
    TestwithPred                -- Test event file with the prediction labels. Generated from SVM
    poslist                     -- The position list of reference.
    Translate the event file to vector. the mod bases were set to 1, unmod bases to 0, and indel bases to nan
    """    
    readnames=np.unique(TestwithPred['read_name']).tolist()
    matrix={}
    b=0
    print(f'There are {len(readnames)} reads mapped on reference')
    for i in readnames:
        b+=1
        if not (b%1000):
            print(b)    #显示进度
        dfi=TestwithPred[TestwithPred['read_name']==i]
        unmodpos=dfi[dfi['Predict']==1]['position'].tolist()
        modpos=dfi[dfi['Predict']==-1]['position'].tolist()
        matrix[i]=[1 if (k in modpos) else (0 if (k in unmodpos) else np.nan) for k in poslist] 
    return pd.DataFrame(matrix).T 

def get_index(lst=None, item=''):
    return [index for (index,value) in enumerate(lst) if value == item]

def Normal_SHAPE(inlist,windowsize=200,step=5): #采用icSHAPE-pipe的策略，在Windows内，top5%为1，最低5%为0
    """
    inlist                      -- A list of rate scores
    windowsize                  -- The sliding window size.
    step                        -- The sliding steps.
    Normalize the SHAPE reactivity scores 
    """    
    shape=[[] for i in range(len(inlist))]
    finalshape=[0 for i in range(len(inlist))]
    start=0
    while start+windowsize< len(inlist):
        sublist=np.array(inlist[start:start+windowsize])
        s95=stats.scoreatpercentile(sublist[~np.isnan(sublist)], 95)
        s5=stats.scoreatpercentile(sublist[~np.isnan(sublist)], 5)
        subshape=[1 if i>s95 else (0 if i<s5 else (i-s5)/(s95-s5)) for i in sublist]
        for i in range(windowsize):
            shape[start+i].append(subshape[i])
        start+=step
    #处理最后一个window
    endsublist=np.array(inlist[start:])
    s95=stats.scoreatpercentile(endsublist[~np.isnan(endsublist)], 95)
    s5=stats.scoreatpercentile(endsublist[~np.isnan(endsublist)], 5)
    subshape=[1 if i>s95 else (0 if i<s5 else (i-s5)/(s95-s5)) for i in endsublist]
    for i in range(len(endsublist)):
        shape[start+i].append(subshape[i])
    #对每个位置的值求平均
    for i in range(len(inlist)):
        finalshape[i]=mean(shape[i])
    return finalshape


def plotScore(shapescore,poslist,outfile):
    """
    shapescore                     -- A list of SHAPE reactivity score
    poslist                        -- A list of the positions.
    outfile                        -- The output figure.
    Plot the SHAPE reactivity scores 
    """    
    plt.figure(figsize=(10,4))
    plt.bar(poslist, shapescore, color =color_SHAPE(shapescore))
    plt.ylabel("Reactivity scores",fontsize=15) 
    plt.xlabel("Position",fontsize=15)
    plt.xlim(0,len(poslist))
    plt.savefig(outfile)
    plt.close()

def plotPCA(pcadata,clust_ids,samplelabel,cluser_ratio,outfile):
    """
    pcadata                       --The PCA dataframe in which include two columns(PC1 and PC2)
    clust_ids                     --A list of cluster labels of each reads. e.g [0,0,1,2,0,1,2]
    cluser_ratio                  -- The reads ratio of each cluster.
    outfile                       -- The output figure.
    Plot the PCA clustering result 
    """    
    cluster_color = { 0: 'red', 1: '#FF7F00', 2: '#FF00FF',3:'#0000FF',4:'yellow',5:'black'}
    fig = plt.figure(figsize = (6,6))
    ax = fig.add_subplot(111)
    plt.scatter(pcadata['PC1'], pcadata['PC2'], c = [cluster_color[x]  if x in samplelabel else 'gray' for x in clust_ids], s = 20, alpha = 0.50 )
    ax.set_xlabel('PC 1', fontsize = 15)
    ax.set_ylabel('PC 2', fontsize = 15)
    N=len(cluser_ratio)
    for i in range(N):
        ax.text(0.05,0.95-0.05*i,'Cluter{}: {:.2%}'.format(i,cluser_ratio[i]),c=cluster_color[i], fontsize = 12, transform=ax.transAxes)
    ax.set_aspect('equal', 'box')
    plt.axis('equal')
    plt.savefig(outfile)
    plt.close()

def color_SHAPE(shape_list, cutoff=[0.3, 0.5, 0.7]):
    """
    shape_list              -- A list of SHAPE scores
    cutoff                  -- Cutoff of SHAPE color boundaries.
    
    Transform SHAPE values to color blocks
    """
    color_blocks=[]
    for value in shape_list:
        if value == 'nan':
            color_blocks.append('lightgray')
        else:
            shape = float(value)
            if shape < cutoff[0]:
                color_blocks.append('black')
            elif shape < cutoff[1]:
                color_blocks.append('blue')
            elif shape < cutoff[2]:
                color_blocks.append('orange')
            else:
                color_blocks.append('red')
    return color_blocks

def CalculateSHAPEFromBitvector_nan(Bitvector):
    Percentlist=[]
    for i in Bitvector.columns:
        s=Bitvector[i]
        if len(s[~np.isnan(s)]) > 0:
            percent=np.nanmean(s) # ignore the None positions and calculate the mean as the rate
            Percentlist.append(float("%.3g"%percent))
        else:
            Percentlist.append(np.nan)
        
    shapescore=Normal_SHAPE(Percentlist)
    return shapescore

def args():
    """
    Calculate and plot the reactivity scores according to different methods
    """ 
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--Mod_Profile", type=str, required=True,      help='The Mod_Profile from SVM.py. Each bases are labeled with the modification information.')
    parser.add_argument("-o", "--output",      type=str, required=True,      help='The output reactivity score file.')
    parser.add_argument("-m", "--method",      type=str, default='mean',     help='The reactivity score calculating methods. mean: calculate the mean score as one single structure; heter: calculate the alternative comformations separately. The default is (mean).')
    parser.add_argument("-l", "--length",   type=int, required=True,      help='The intermidiate length.')
    parser.add_argument("-b", "--bias",    type=int, required=True, help='The intermidiate length error tolerance, reflecting the resolution of the intermidiates. Reads that 3 terminal mapping position ranged in end~end+bias were grouped as one intermidiate.')


    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = args()

    if not os.path.exists(args.Mod_Profile):
        print(f'Mod_Profile file not found: {args.Mod_Profile}')
        exit()
    assert args.method in ('mean', 'heter'), "method should be one of mean/heter"

    #read in the Mod_Profile
    infile=pd.read_csv(args.Mod_Profile,sep='\t',index_col=0)

    #get the intermidiate length
    poslist=[i for i in range(1,args.length+args.bias+1)]
    
    #calculate the scores according to the methods
    if args.method=='mean':
        rates=[]
        for i in poslist:
            temp=infile[infile['position']==i]
            if temp.shape[0]>=1:
                pred=temp['Predict']
                rates.append(np.sum(pred==-1)/temp.shape[0])
            else:
                rates.append(np.nan)
        shapescore=Normal_SHAPE(rates)
        #Plot the PCA and Clustering results into PDF
        plotScore(shapescore,poslist,args.output+'.pdf')
        #Write the scores to outfile
        shapescoredit=pd.DataFrame(shapescore,index=poslist,columns=['score'])
        #shapescoredit.fillna('Null').to_csv(args.output+'.score',sep='\t')
        shapescoredit.fillna('NaN').to_csv(args.output+'.score',sep='\t')
        print("The mean Reactivity scores were calculated!")


    if args.method=='heter':
        
        #transfer or read the profile to Bitvector
        if not os.path.exists(args.output+'.vect'):
            vect=transPredictedEventToMatrix(infile,poslist)
            vect.to_csv(args.output+'.vect',sep='\t',header=True, index=True)
            print("The Bitvector was converted!\n")
        else:
            vect=pd.read_csv(args.output+'.vect',sep='\t',index_col=0)
            print("The Bitvector was found!\n")

        # PCA
        pca = PCA(n_components=2)
        temp_vect=vect.dropna(axis=1,how='all') #exclude the columns that are all None
        pca_out=pd.DataFrame(pca.fit_transform(temp_vect.fillna(0)),columns = ['PC1', 'PC2'])
        
        plt.scatter(pca_out['PC1'], pca_out['PC2'])
        plt.savefig(args.output+'.PCA.pdf')
        plt.close()

        # Check the PCA.PDF file and assign the cluster number and method
        N=input("Please check the PCA.pdf file and then assign the Number of Clusers:  ")
        print(f"Number of Clusers is assigned:{N}!\n")

        M=input("Please choose the Clustering method, G for GaussianMixture and D for DBSCAN:  ")
        #assert M in ('G', 'D'), "method should be one of G/D"
        print(f"The cluser method is assigned:{M}!\n")

        if M=='G':
            clustering =GaussianMixture(n_components = int(N),covariance_type='full',random_state=1).fit(pca_out)
            cluster_ids=clustering.predict(pca_out).tolist()
            lists=Counter(TPP_clust_ids).most_common(int(N))
            samplelabel={int(k):v for k, v in lists}.keys() 

        if M=='D':
            clustering =DBSCAN(eps=0.1, min_samples=int(N)).fit(pca_out)
            cluster_ids =clustering.labels_.tolist()
            lists=Counter(TPPclustering.labels_).most_common(int(N))
            samplelabel={int(k):v for k, v in lists}.keys() 

        #Calculate the scores of each cluster
        cluser_ratio=[]
        cluser_vect=[]
        cluser_score=[]
        cluser_name=[]
        for i in range(int(N)):
            cluser_ratio.append(cluster_ids.count(i)/len(cluster_ids))
            cluser_vect.append(vect.iloc[get_index(cluster_ids,i),:])
            cluser_score.append(CalculateSHAPEFromBitvector_nan(cluser_vect[i]))
            cluser_name.append('Cluster'+str(i+1))
        print("Scores calculate successfully!")
        #Plot the PCA and Clustering results into PDF
        plotPCA(pca_out,cluster_ids,samplelabel,cluser_ratio,args.output+'.PCA.pdf')
        #Write the scores to outfile
        shapescoredit=pd.DataFrame(cluser_score,index=cluser_name,columns=poslist).T
        shapescoredit.fillna('Null').to_csv(args.output+'.scors',sep='\t')
        print("The heterogeneous reactivity scores were calculated!")

