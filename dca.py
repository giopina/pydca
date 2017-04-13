# Direct Coupling Analysis (DCA)
#
# function dca(inputfile , outputfile)
# 
# INPUTS: 
#   inputfile  - file containing the FASTA alignment
#   outputfile - file for dca results. The file is composed by N(N-1)/2 
#                (N = length of the sequences) rows and 4 columns: 
#                residue i (column 1), residue j (column 2),
#                MI(i,j) (Mutual Information between i and j), and 
#                DI(i,j) (Direct Information between i and j).
#                Note: all insert columns are removed from the alignment.
#
# SOME RELEVANT VARIABLES:
#   N        number of residues in each sequence (no insert)
#   M        number of sequences in the alignment
#   Meff     effective number of sequences after reweighting
#   q        equal to 21 (20 aminoacids + 1 gap)
#   align    M x N matrix containing the alignmnent
#   Pij_true N x N x q x q matrix containing the reweigthed frequency
#            counts.
#   Pij      N x N x q x q matrix containing the reweighted frequency 
#            counts with pseudo counts.
#   C        N(q-1) x N(q-1) matrix containing the covariance matrix.
#
#
# Copyright for this implementation: 
#             2011/12 - Andrea Pagnani and Martin Weigt
#                       andrea.pagnani@gmail.com 
#                       martin.weigt@upmc.fr
# 
# Permission is granted for anyone to copy, use, or modify this
# software and accompanying documents for any uncommercial
# purposes, provided this copyright notice is retained, and note is
# made of any changes that have been made. This software and
# documents are distributed without any warranty, express or
# implied. All use is entirely at the user's own risk.
#
# Any publication resulting from applications of DCA should cite:
#
#     F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander, 
#     R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
#     analysis of residue co-evolution captures native contacts across 
#     many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np
import support_functions as sf

class DCA:
    """Class DCA:
    direct coupling analysis
    I don't know yet how it will work"""
    
    def __init__(self,inputfile,pseudocount_weight=0.5,theta=0.1):
        """Constructor of the class DCA"""
        self.pseudocount_weight=pseudocount_weight
        self.theta=theta

        fasta_list=sf.FASTA_parser(inputfile,check_aminoacid=True)
        self.alignment=sf.Alignment(fasta_list)
        self.N=self.alignment.N
        self.M=self.alignment.M
        
        self.Compute_True_Frequencies()

    def Compute_True_Frequencies(self):
        """Computes reweighted frequency counts"""
        ### TODO: this still has to be checked and tested
        from scipy.spatial.distance import pdist
        from scipy.spatial.distance import squareform
        q=self.alignment.q
        W = np.ones(self.M)
        align=self.alignment.Z
        if self.theta > 0.0 :
            #   W = (1./(1+sum(squareform(pdist(align,'hamming')<theta))));
            cacca=(pdist(align,metric='hamming')<self.theta)
            print(cacca.shape)
            W= (1./(1+np.sum(squareform(cacca),axis=0)))
            print(W.shape)
        self.Meff=np.sum(W)

        self.Pij_true = np.zeros((self.N,self.N,q,q))
        self.Pi_true = np.zeros((self.N,q))
        import time
####### this is the original way of doing it but mine is faster
#        t0=time.time()
#        for j in range(self.M):
#            for i in range(self.N):
#                self.Pi_true[i,align[j,i]] = self.Pi_true[i,align[j,i]] + W[j]
#        self.Pi_true = self.Pi_true/self.Meff
#        print(self.Pi_true,time.time()-t0)
############################################
########## this way is x10 faster
        self.Pi_true=np.zeros((self.N,q))
        #        for i in range(self.N):
        t0=time.time()
        for a in range(q):
            self.Pi_true[:,a]=np.sum(((align==a)*W[:,np.newaxis]),axis=0)
            #            Pi_true[:,align[j,:]]+=W[j]
        self.Pi_true/=self.Meff
##############################
#        print(Pi_true,time.time()-t0)
#        print(np.allclose(self.Pi_true,Pi_true))

#        t0=time.time()
###### this is the original way but mine is faster (at least for the alignment sizes I'm dealing with (700x500)
#        for l in range(self.M):
#            for i in range(self.N-1):
#                for j in range(i+1,self.N):
#                    self.Pij_true[i,j,align[l,i],align[l,j]] = self.Pij_true[i,j,align[l,i],align[l,j]] + W[l]
#                    # self.Pij_true[j,i,align[l,j],align[l,i]] = self.Pij_true[i,j,align[l,i],align[l,j]] ### no sense in doing this here. I can simmetrize later
#        self.Pij_true+=np.transpose(self.Pij_true,axes=(1,0,3,2))
#        self.Pij_true = self.Pij_true/self.Meff
#
#        ### I think the following triple loop can be simplyfied a lot...
#        scra = np.eye(q);
#        for i in range(self.N):
#            for alpha in range(q):
#                for beta in range(q):
#                    self.Pij_true[i,i,alpha,beta] = self.Pi_true[i,alpha] * scra[alpha,beta]
#        ###
#        print(time.time()-t0)
#####################################33
        self.Pij_true=np.zeros((self.N,self.N,q,q))
        #t0=time.time()
        for a in range(q):
            for b in range(q):
                self.Pij_true[:,:,a,b]+=np.tensordot((align==a)*W[:,np.newaxis],(align==b),axes=(0,0))
        self.Pij_true = self.Pij_true/self.Meff
        #print(time.time()-t0)
        #print(np.allclose(self.Pij_true,Pij_true))
            
        

#    return Pij_true,Pi_true,Meff

###
### Here I'm just copying the headers of the functions defined in dca.m
###

#def return_alignment(inputfile):
#    """Reads alignment from inputfile, removes insters and converts into numbers"""
#    return N,M,q,Z
### now this is done via the class Alignment (see support_functions.py)

    
#def Compute_True_Frequencies(align,M,N,q,theta):
#    """Computes reweighted frequency counts"""
#    return Pij_true,Pi_true,Meff


def Compute_Results(Pij,Pi,Pij_true,Pi_true,invC,N,q,fp):
    """Computes and prints the mutual and direct informations"""

def with_pc(Pij_true,Pi_true,pseudocount_weight,N,q):
    """Adds pseudocounts"""
    return Pij,Pi

def Compute_C(Pij,Pi,N,q):
    """Computes correlation matrix"""
    return C

def mapkey(i,alpha,q):
    # not sure what this is doing
    return A

def calculate_mi(i,j,P2,P1,q):
    """Computes mutual information between columns i and j"""
    return M,s1,s2

def ReturnW(C,i,j,q):
    """Extracts coupling matrix for columns i and j"""
    return W

def bp_link(i,j,W,P1,q):
    """Computes direct information"""
    return DI

def compute_mu(i,j,W,P1,q):
    # not sure what this is doing
    return mu1,mu2

def compute_di(i,j,W, mu1,mu2, Pia):
    """computes direct information"""
    return DI


