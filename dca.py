# Direct Coupling Analysis (DCA)
#
# 
# SOME RELEVANT VARIABLES:
#   N          number of residues in each sequence (no insert)
#   M          number of sequences in the alignment
#   Meff       effective number of sequences after reweighting
#   q          equal to 21 (20 aminoacids + 1 gap)
#   alignment  M x N matrix containing the alignmnent
#   Pij_true   N x N x q x q matrix containing the reweigthed frequency
#              counts.
#   Pij        N x N x q x q matrix containing the reweighted frequency 
#              counts with pseudo counts.
#   invC          N x (q-1) x N x (q-1) inverse of the covariance matrix.
#
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
    direct coupling analysis"""
    
    #def __init__(self,inputfile,pseudocount_weight=0.5,theta=0.1,get_MI=False,get_DI=True):
    def __init__(self,alignment,pseudocount_weight=0.5,theta=0.1,get_MI=False,get_DI=True):
        """Constructor of the class DCA"""
        self.pseudocount_weight=pseudocount_weight # relative weight of pseudo count
        self.theta=theta # threshold for sequence id in reweighting

        #fasta_list=sf.FASTA_parser(inputfile,check_aminoacid=True)
        #self.alignment=sf.Alignment(fasta_list)
        self.alignment=alignment
        self.N=self.alignment.N
        self.M=self.alignment.M
        self.q=self.alignment.q
        print("compute true frequencies...")
        self.__comp_true_freq()
        if get_MI:
            print("compute MI")
            self.__comp_MI()
        print("add pseudocounts")
        self.__with_pc()
        print("compute C")
        self.__comp_C()
        if get_DI:
            print("compute DI")
            self.__comp_DI()
        print("Done!")

    def __comp_true_freq(self):
        """Computes reweighted frequency counts"""
        from scipy.spatial.distance import pdist
        from scipy.spatial.distance import squareform
        W = np.ones(self.M)
        align=self.alignment.Z
        if self.theta > 0.0 :
            cacca=(pdist(align,metric='hamming')<self.theta)
            W= (1./(1+np.sum(squareform(cacca),axis=0)))
        self.Meff=np.sum(W)

        self.__Pij = np.zeros((self.N,self.N,self.q,self.q))
        self.__Pi = np.zeros((self.N,self.q))

        for a in range(self.q):
            self.__Pi[:,a]=np.sum(((align==a)*W[:,np.newaxis]),axis=0)
        self.__Pi/=self.Meff

        for a in range(self.q):
            for b in range(self.q):
                self.__Pij[:,:,a,b]+=np.tensordot((align==a)*W[:,np.newaxis],(align==b),axes=(0,0))
        self.__Pij = self.__Pij/self.Meff
            
    def __with_pc(self):
        """Adds pseudocounts"""
        self.__Pi = (1.-self.pseudocount_weight)*self.__Pi +\
                  self.pseudocount_weight/self.q*np.ones((self.N,self.q))
        ### TODO: do we need to store both __Pij_true and __Pij?
        Pij_diag=self.__Pij[range(self.N),range(self.N),:,:]
        self.__Pij = (1.-self.pseudocount_weight)*self.__Pij +\
                   self.pseudocount_weight/self.q/self.q*np.ones((self.N,self.N,self.q,self.q))
        scra = np.eye(self.q)
        for i in range(self.N):
            self.__Pij[i,i,:,:] = (1.-self.pseudocount_weight)*Pij_diag[i,:,:] +\
                            self.pseudocount_weight/self.q*scra

    def __comp_C(self):
        """Computes correlation matrix and its inverse"""
        ### Remember remember... I'm excluding A,B = q (see PNAS SI, pg 2, column 2)
        C=np.transpose(\
                            self.__Pij[:,:,:-1,:-1] -\
                            self.__Pi[:,np.newaxis,:-1,np.newaxis]*\
                            self.__Pi[np.newaxis,:,np.newaxis,:-1],\
                            axes=(0,2,1,3))
        del self.__Pij
        #del self.__Pi
        ### NB: the order of the indexes in C is different from __Pij, this is needed for tensorinv. TODO: think if it's better to use the same order of indexes for every array
        from numpy.linalg import tensorinv
        self.__invC=tensorinv(C)

    def __comp_MI(self):
        """Computes the mutual information"""
        ### TODO: there should be a smarter way of storing and accessing these symmetric matrices using only half the space
        self.mutual_information=np.zeros((self.N,self.N))
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                # mutual information
                self.mutual_information[i,j] = self.__calc_mi(i,j)
        self.mutual_information+=self.mutual_information.T

            
    def __calc_mi(self,i,j):
        """Computes mutual information between columns i and j"""
        ### Here apparently I'm using Pij_true
        ### Apparently, also I do not need this useless Mutual information...
        ### Even more apparently, two of the three output of this function are not even used in the matlab code (s1,s2, aka si_true, sj_true)....
        M = 0.
        ### TODO: this loop can probably be rewritten in a smart "numpy" way
        for alpha in range(self.q):
            for beta in range(self.q):
                if self.__Pij[i,j,alpha,beta]>0:
                    M = M + self.__Pij[i,j,alpha, beta]*np.log(self.__Pij[i,j, alpha, beta] / self.__Pi[i,alpha]/self.__Pi[j,beta])
        return M
    
    def Print_Results(self,filename):
        fh=open(filename,'w')
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                fh.write('%d %d '%(i, j))                
                fh.write('%g '%self.mutual_information)
                fh.write('%g '%self.direct_information)
                fh.write('\n')
        fh.close()
          
    def __comp_DI(self):
        """Computes Direct Information"""
        ### TODO: implement other DCAmethods
        ### TODO: there should be a smarter way of storing and accessing these symmetric matrices using only half the space
        self.direct_information=np.zeros((self.N,self.N))
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                # direct information from mean-field
                self.direct_information[i,j] = self.__bp_link(i,j)
        self.direct_information+=self.direct_information.T
        del __invC
        del self.__Pi
        
    def __bp_link(self,i,j):
        """Computes direct information"""
        W_mf=np.ones((self.q,self.q))
        W_mf[:-1,:-1]= np.exp( -self.__invC[i,:,j,:] )
        mu1, mu2 = self.__comp_mu(i,j,W_mf);
        DI = self.__comp_di(i,j,W_mf, mu1,mu2);
        return DI
    
    def __comp_mu(self,i,j,W):
        ### not sure what this is doing
        epsilon=1e-4
        diff =1.0
        mu1 = np.ones((1,self.q))/self.q
        mu2 = np.ones((1,self.q))/self.q
        pi = self.__Pi[i,:]
        pj = self.__Pi[j,:]

        while ( diff > epsilon ):
            ### TODO: add a counter and a maxiter parameter?
            scra1 = np.dot(mu2, W.T)
            scra2 = np.dot(mu1, W)
            new1 = pi/scra1
            new1 /=np.sum(new1)
            new2 = pj/scra2
            new2 /= np.sum(new2)

            diff = max( (np.abs( new1-mu1 )).max(), (np.abs( new2-mu2 )).max() )
            mu1 = new1
            mu2 = new2
        return mu1,mu2

    def __comp_di(self,i,j,W, mu1,mu2):
        """computes direct information"""
        tiny = 1.0e-100
        Pdir = W*np.dot(mu1.T,mu2)
        Pdir = Pdir / np.sum(Pdir)
        Pfac = self.__Pi[i,:][:,np.newaxis]*self.__Pi[j,:][np.newaxis,:]
        ### TODO why trace? Shouldn't it be the sum over all elements?
        DI = np.trace(\
                      np.dot(Pdir.T , np.log((Pdir+tiny)/(Pfac+tiny)) ) \
        )
        return DI

    def get_ordered_di(self,k_pairs=None,offset=4,return_di=False):
        """Sort the pairs by their direct information"""
        if k_pairs==None:
            k_pairs=self.N*2
        faraway=np.triu_indices(self.N,k=offset)
        self.di_order=(np.array(faraway).T[np.argsort(self.direct_information[faraway])])[::-1]
        if return_di:
            return self.direct_information[[self.di_order[:k_pairs,0],self.di_order[:k_pairs,1]]] ### TODO: this is not creating a copy. Be careful


def print_contacts(dca_obj,n_di=None,colore='b',lower_half=False):
    """Prints the contact maps derived from a DCA object"""
        #plt.figure(figsize=(8,8))
        #n_di=1000
    if n_di==None:
        n_di=dca_obj.N*2
    if lower_half:
        plt.scatter(dca_obj.di_order[:n_di,1],dca_obj.di_order[:n_di,0],\
                    marker='s',s=3,color=colore)
    else:
        plt.scatter(dca_obj.di_order[:n_di,0],dca_obj.di_order[:n_di,1],\
                    marker='s',s=3,color=colore)
    plt.plot(range(dca_obj.N),color='black')
    plt.ylim(0,dca_obj.N)
    plt.xlim(0,dca_obj.N)


def compute_dca(inputfile,pseudocount_weight=0.5,theta=0.1):
    """Perform mfDCA starting from a FASTA input file. Returns a DCA object"""
    ### TODO: add "filter" argument to filter sequences with too many gaps. add "method" arguments to use different DCA implementations
    fasta_list=sf.FASTA_parser(inputfile,check_aminoacid=True)
    alignment=sf.Alignment(fasta_list)
    dca_obj=dca.DCA(alignment,pseudocount_weight=pseudocount_weight,theta=theta,get_MI=False,get_DI=True)
    dca_obj.get_ordered_DI()
    return dca_obj
