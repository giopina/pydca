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
        self.pseudocount_weight=pseudocount_weight # relative weight of pseudo count
        self.theta=theta # threshold for sequence id in reweighting

        fasta_list=sf.FASTA_parser(inputfile,check_aminoacid=True)
        self.alignment=sf.Alignment(fasta_list)
        self.N=self.alignment.N
        self.M=self.alignment.M
        self.q=self.alignment.q
        print("compute true frequencies...")
        self.Compute_True_Frequencies()
        print("add pseudocounts")
        self.with_pc()
        print("compute C")
        self.Compute_C()
        print("compute results")
        self.Compute_Results('prova-DI.dat')
        print("Done!")

    def Compute_True_Frequencies(self):
        """Computes reweighted frequency counts"""
        from scipy.spatial.distance import pdist
        from scipy.spatial.distance import squareform
        W = np.ones(self.M)
        align=self.alignment.Z
        if self.theta > 0.0 :
            #   W = (1./(1+sum(squareform(pdist(align,'hamming')<theta))));
            cacca=(pdist(align,metric='hamming')<self.theta)
            print(cacca.shape)
            W= (1./(1+np.sum(squareform(cacca),axis=0)))
            print(W.shape)
        self.Meff=np.sum(W)

        self.Pij_true = np.zeros((self.N,self.N,self.q,self.q))
        self.Pi_true = np.zeros((self.N,self.q))
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
        self.Pi_true=np.zeros((self.N,self.q))
        #        for i in range(self.N):
#        t0=time.time()
        for a in range(self.q):
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
        self.Pij_true=np.zeros((self.N,self.N,self.q,self.q))
        #t0=time.time()
        for a in range(self.q):
            for b in range(self.q):
                self.Pij_true[:,:,a,b]+=np.tensordot((align==a)*W[:,np.newaxis],(align==b),axes=(0,0))
        self.Pij_true = self.Pij_true/self.Meff
        #print(time.time()-t0)
        #print(np.allclose(self.Pij_true,Pij_true))
            
    def with_pc(self):
        """Adds pseudocounts"""
        ### TODO: do we need to store both Pij_true and Pij?
        ### TODO: not tested
        self.Pij = (1.-self.pseudocount_weight)*self.Pij_true +\
                   self.pseudocount_weight/self.q/self.q*np.ones((self.N,self.N,self.q,self.q))
        self.Pi = (1.-self.pseudocount_weight)*self.Pi_true +\
                  self.pseudocount_weight/self.q*np.ones((self.N,self.q))
        Pij=np.array(self.Pij)
        scra = np.eye(self.q)
        for i in range(self.N):
            self.Pij[i,i,:,:] = (1.-self.pseudocount_weight)*self.Pij_true[i,i,:,:] +\
                            self.pseudocount_weight/self.q*scra

    def Compute_C(self):
        """Computes correlation matrix"""
        ### TODO: not tested
        ### Remember remember... I'm excluding A,B = q (see PNAS SI, pg 2, column 2)

        ### che e' sta roba?? suppongo si possa fare meglio con qualche reshape e operazioni fra array...
        #self.C=np.zeros((self.N*(self.q-1),self.N*(self.q-1)))
        #for i in range(self.N):
        #    for j in range(self.N):
        #        for alpha in range(self.q-1):
        #            for beta in range(self.q-1):
        #                self.C[self.mapkey(i,alpha),self.mapkey(j,beta)] = self.Pij[i,j,alpha,beta] - self.Pi[i,alpha]*Pi[j,beta]
        ##############################
        #    #def mapkey(self,i,alpha):
        #    #    A = (self.q-1)*i+alpha
        #    #    return A
        ##############################

        ### mi ghe provo
#        cacca=self.Pi[:,np.newaxis,:,np.newaxis]*self.Pi[np.newaxis,:,np.newaxis,:]
#        print(cacca.shape,self.Pij.shape)
        self.C=np.transpose(\
                            self.Pij[:,:,:-1,:-1] -\
                            self.Pi[:,np.newaxis,:-1,np.newaxis]*\
                            self.Pi[np.newaxis,:,np.newaxis,:-1],\
                            axes=(0,2,1,3))
        ### NB: the order of the indexes in C is different from Pij, this is needed for tensorinv. TODO: think if it's better to use the same order of indexes for every array
        from numpy.linalg import tensorinv
        self.invC=tensorinv(self.C)


    def Compute_Results(self,filename):
        """Computes and prints the mutual and direct informations"""
        fh=open(filename,'w')
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                # mutual information
                MI_true = self.calculate_mi(i,j);
                
                # direct information from mean-field
                W_mf=np.ones((self.q,self.q))
                W_mf[:-1,:-1]= np.exp( -self.invC[i,:,j,:] ) #self.ReturnW(i,j);                
                DI_mf_pc = self.bp_link(i,j,W_mf);
                fh.write('%d %d %g %g\n'%(i, j, MI_true, DI_mf_pc))
        fh.close()
            
    def calculate_mi(self,i,j):
        """Computes mutual information between columns i and j"""
        ### Here apparently I'm using Pij_true
        ### Apparently, also I do not need this useless Mutual information...
        ### Even more apparently, two of the three output of this function are not even used in the matlab code (s1,s2, aka si_true, sj_true)....
        M = 0.
        for alpha in range(self.q):
            for beta in range(self.q):
                if self.Pij_true[i,j,alpha,beta]>0:
                    M = M + self.Pij_true[i,j,alpha, beta]*np.log(self.Pij_true[i,j, alpha, beta] / self.Pi_true[i,alpha]/self.Pi_true[j,beta])
                        
            #s1=0.
            #s2=0.
            #for alpha in range(q):
            #    if( self.Pi_true[i,alpha]>0 ):
            #        s1 = s1 - self.Pi_true[i,alpha] * np.log(self.Pi_true[i,alpha])
            #    if( self.Pi_true[j,alpha]>0 ):
            #        s2 = s2 - self.Pi_true[j,alpha] * np.log(self.Pi_true[j,alpha])
        return M#,s1,s2

    def bp_link(self,i,j,W_mf):
        """Computes direct information"""
        mu1, mu2 = self.compute_mu(i,j,W_mf);
        DI = self.compute_di(i,j,W_mf, mu1,mu2);
        return DI
    
    def compute_mu(self,i,j,W):
        ### not sure what this is doing
        epsilon=1e-4
        diff =1.0
        mu1 = np.ones((1,self.q))/self.q
        mu2 = np.ones((1,self.q))/self.q
        pi = self.Pi[i,:]
        pj = self.Pi[j,:]
        #print(mu1.shape,W.shape,pi.shape)
        while ( diff > epsilon ):
            ### TODO: add a counter and a maxiter parameter!
            scra1 = np.dot(mu2, W.T)
            scra2 = np.dot(mu1, W)
            new1 = pi/scra1
            new1 /=np.sum(new1)
            new2 = pj/scra2
            new2 /= np.sum(new2)

            ### TODO: check this diff. I'm not sure I'm doing the right max
            diff = max( (np.abs( new1-mu1 )).max(), (np.abs( new2-mu2 )).max() )

            mu1 = new1
            mu2 = new2
        #print(mu1)
        return mu1,mu2

    def compute_di(self,i,j,W, mu1,mu2):
        """computes direct information"""
        tiny = 1.0e-100
        #Pdir = W*(mu1[np.newaxis,:]*mu2[:,np.newaxis])
        Pdir = W*np.dot(mu1.T,mu2)
        Pdir = Pdir / np.sum(Pdir)
        #print(Pdir.shape)
        Pfac = self.Pi[i,:][np.newaxis,:]*self.Pi[j,:][:,np.newaxis]
        #print(Pfac.shape)
        ### why trace? Shouldn't it be the sum over all elements?
        DI = np.trace(\
                      np.dot(Pdir.T , np.log((Pdir+tiny)/(Pfac+tiny)) ) \
        )
#        DI = np.trace(\
#                    Pdir[np.newaxis,:]*\
#                      np.log( (Pdir+tiny)/(Pfac+tiny) )\
#            )
        
        return DI

        
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

#def with_pc(Pij_true,Pi_true,pseudocount_weight,N,q):
#    """Adds pseudocounts"""
#    return Pij,Pi

#def Compute_Results(Pij,Pi,Pij_true,Pi_true,invC,N,q,fp):
#    """Computes and prints the mutual and direct informations"""

#def Compute_C(Pij,Pi,N,q):
#    """Computes correlation matrix"""
#    return C

#def mapkey(i,alpha,q):
#    # not sure what this is doing
#    return A

#def calculate_mi(i,j,P2,P1,q):
#    """Computes mutual information between columns i and j"""
#    return M,s1,s2

#def ReturnW(C,i,j,q):
#    """Extracts coupling matrix for columns i and j"""
#    return W

#def bp_link(i,j,W,P1,q):
#    """Computes direct information"""
#    return DI

#def compute_mu(i,j,W,P1,q):
#    # not sure what this is doing
#    return mu1,mu2

#def compute_di(i,j,W, mu1,mu2, Pia):
#    """computes direct information"""
#    return DI


