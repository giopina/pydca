import numpy as np
import matplotlib.pyplot as plt
import pydca.msa as msa
import pydca.dca as dca

def plot_gap_in_seq(alin):
    gap_in_seq=np.sum(alin.Z==0,axis=1)

    plt.figure(figsize=(13,4))
    plt.subplot2grid((1,2),(0,0))
    plt.plot(np.sort(gap_in_seq),marker='',ls='-')
    plt.ylabel('# of gaps in sequence')
    plt.xlabel('sequence index')
    
    plt.subplot2grid((1,2),(0,1))
    cacca=plt.hist(gap_in_seq,bins=np.arange(0,alin.N,20),histtype='bar',edgecolor='magenta',color='0.9')
    plt.yscale('log')
    cacca=plt.xticks(np.arange(0,alin.N,20))
    plt.ylabel('N. of sequences')
    plt.xlabel('N. of gaps')


def plot_gap_per_pos(alin):
    gap_per_pos=np.sum(alin.Z==0,axis=0)
    
    plt.figure(figsize=(15,3))
    plt.bar(np.arange(alin.N),gap_per_pos,width=1.0,color='green')
    plt.xlim(-1,alin.N)
    plt.ylabel('N. of gaps per position')
    plt.xlabel('aa index')
    
    plt.figure(figsize=(15,3))
    plt.bar(np.arange(alin.N),gap_per_pos/alin.M,width=1.0,color='red')
    plt.xlim(-1,alin.N)
    plt.ylim(0,1)
    plt.ylabel('Freq. of gaps per position')
    plt.xlabel('aa index')
    plt.grid()

def plot_gap_len(alin):
    import re
    plt.figure(figsize=(12,4))
    stringhe_di_gap=[re.sub("[A-Z]",'/',s).split('/') for s in alin.stripped_seqs]
    gap_len=[len(s) for seq in stringhe_di_gap for s in seq if len(s)>0 ]
    plt.subplot2grid((1,2),(0,0))
    cacca=plt.hist(gap_len,bins=np.arange(0,alin.N,10),histtype='step')
    plt.yscale('log')
    
    n,b=np.histogram(gap_len,bins=np.arange(0,alin.N,10))
    plt.subplot2grid((1,2),(0,1))
    plt.plot(b[:-1],n/alin.M,ls='',marker='x')
    plt.yscale('log')

def insert_gaps(alin,gap_len,gap_num,seed=12345,save_seq=None):
    """
    This function inserts stretches of gaps in a given alignment. 
    gap_len is the length of the gap stretches you want to insert.
    gap_num is how many of these stretches you want to insert (so total number of gaps is actually gap_len*gap_num):
    It inserts the stretches in random sequences at random positions.
    Returns a new alignment object.
    save_seq should be an index of a sequence you don't want to modify (if any).
    """
    fasta_list=alin.get_dict(stripped=True)
    np.random.seed(12345)
    for k in range(gap_num):
        iseq=np.random.randint(0,high=alin.M)
        if iseq==save_seq:
            continue
        ipos=np.random.randint(0,high=alin.N-gap_len)
        fasta_list[iseq]['sequence']=fasta_list[iseq]['sequence'][:ipos]+\
                                      '-'*gap_len+fasta_list[iseq]['sequence'][ipos+gap_len:]
    new_alignment=msa.Alignment(fasta_list)
    return new_alignment
