import argparse
import os
import pandas as pd
import numpy as np


def select_intermidiates_from_event_single(Eventfile,gene,outdir,start,end,bias): #从指定gene中按照reads的mapping位置分离出reads
    if not os.path.exists(outdir):
        os.system(f'mkdir {outdir}')
    infile=pd.read_csv(Eventfile,sep='\t',index_col=0,header=(0))
    df=infile[infile['contig']==gene]
    print('Open event file finished...')
    
    outfilename=outdir+'/'+gene+'_'+str(end)+'.event'
    outfile=open(outfilename,'a+')
    head=df.columns.tolist()
    outfile.writelines('\t'.join(head)+'\n') #打印表头
    count=0
    for i in df.groupby('read_name'): #遍历每个read，按照长度分类
        pos=i[1]['position'].tolist()
        EndPos=pos[-1]
        StartPos=pos[0]
        if StartPos<start and EndPos>=end and EndPos<=end+bias:
            i[1].to_csv(outfile,sep='\t',header=False, index=False)
            count+=1
    outfile.close()
    return count


def select_intermidiates_from_event_multi(Eventfile,gene,outdir,start,end,bias): 
    if not os.path.exists(outdir):
        os.system(f'mkdir {outdir}')
    infile = pd.read_csv(Eventfile,sep='\t',index_col=0,header=(0)) 
    df=infile[infile['contig']==gene]
    print('Reading Event-file finished...')

    n=int(end/bias) #计算一共要分成多少个群
    outlist=[str(i) for i in range(n+1)]  #输出文件的名称list
    countlist=[0 for i in range(n+1)] #每个分群的reads条数list

    head=df.columns.tolist()
    for i in range(n+1):
        outfilename=outdir+'/'+gene+'_'+str(bias*i)+'.event'
        outlist[i]=open(outfilename,'a+')
        outlist[i].writelines('\t'.join(head)+'\n')#打印表头
    print('Open split files finished...')
    for i in df.groupby('read_name'): #遍历每个read，按照长度分类
        pos=i[1]['position'].tolist()
        EndPos=pos[-1]
        StartPos=pos[0]
        if StartPos<start and EndPos>=start and EndPos<=end+bias:
            x=int(EndPos/bias)
            i[1].to_csv(outlist[x],sep='\t',header=False, index=False) 
            countlist[x]+=1
    for i in range(n+1):
        outlist[i].close()
    return countlist


def args():
    """
    select RNA intermidiates based on the mapped read's start and end positions on reference

    """ 
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--event",   type=str, required=True,      help='The event file from dataprocessing.py.')
    parser.add_argument("-g", "--gene", type=str, required=True,      help='The specific gene name in the reference.')
    parser.add_argument("-d", "--output_dir",  type=str, required=True,     help='The output file path where to store the selected intermidiates.')
    parser.add_argument("-s", "--start",   type=int, required=True,      help='The maximal 5 terminal  position of a intermidiate. Reads that 5 terminal mapped beyond the start site were excluded.')
    parser.add_argument("-e", "--end",    type=int, required=True,  help='The 3 terminal position of a intermidiates. Reads that 3 terminal mapped within the end site were excluded.')
    parser.add_argument("-b", "--bias",    type=int, required=True, help='The intermidiate length error tolerance, reflecting the resolution of the intermidiates. Reads that 3 terminal mapping position ranged in end~end+bias were grouped as one intermidiate.')
    parser.add_argument("-m", "--multi",  action='store_true', default=False, help='If the parameter was used, the all intermidiates within the end position were selected parallelly.')

    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = args()

    if not os.path.exists(args.event):
        print(f'Event file not found: {args.event}')
        exit()
    if not os.path.exists(args.output_dir):
        os.system(f'mkdir {args.output_dir}')
    
    Eventfile=args.event
    gene=args.gene
    outdir=args.output_dir
    start=args.start
    end=args.end
    bias=args.bias

    if args.multi:
        countlist=select_intermidiates_from_event_multi(Eventfile,gene,outdir,start,end,bias)
        n=len(countlist)
        print (f'The number of each intermidiates {countlist}!')
    else:
        count=select_intermidiates_from_event_single(Eventfile,gene,outdir,start,end,bias)
        print (f'There are {count} reads in the selected intermidiate!')




