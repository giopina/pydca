# Direct Coupling Analysis (DCA) in python
#
# This code implements the mean field DCA method
#
# Written by G. Pinamonti and E. Sarti
# @ FU Berlin and NIH, in 2017
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
#   invC          N x (q-1) x N x (q-1) inverse of the covariance matrix
#
# Parts of this code are based on the Matlab implementation freely
# available for download at http://dca.rice.edu/portal/dca/
#
#
# Permission is granted for anyone to copy, use, or modify this
# software and accompanying documents for any uncommercial
# purposes, provided this copyright notice is retained, and note is
# made of any changes that have been made. This software and
# documents are distributed without any warranty, express or
# implied. All use is entirely at the user's own risk.
#
# Any publication resulting from applications of mfDCA should cite:
#
#     F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander, 
#     R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
#     analysis of residue co-evolution captures native contacts across 
#     many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np
import support_functions as sf
import matplotlib.pyplot as plt # TODO: where should this go? here or in the plot function?

class DCA:
    """Class DCA:
    direct coupling analysis"""
    
    def __init__(self,alignment,pseudocount_weight=0.5,theta=0.1,get_MI=False,get_DI=True,get_CFN=False):
        """Constructor of the class DCA"""
        ### TODO: what if we put a "save_invC" flag? suppose one wants to work with the couplings afterwards...
        
        self.pseudocount_weight=pseudocount_weight # relative weight of pseudo count
        self.theta=theta # threshold for sequence id in reweighting
        
        self.alignment=alignment
        self.N=self.alignment.N
        self.M=self.alignment.M
        self.q=self.alignment.q
        print("compute true frequencies...")
        self.__comp_true_freq()
        if get_MI:
            print("compute mutual information, MI")
            self.__comp_MI()
        print("add pseudocounts")
        # pseudocounts to avoid singularities due to elements equal to zero
        self.__add_pseudocounts()
        print("compute correlation matrix, C")
        self.__comp_C()
        
        if get_DI:
            print("compute direct information, DI")
            self.__comp_DI()
            self.get_ordered_di() ### Let's always do it
        print("Done!")
        if get_CFN:
            print("compute corrected Frobenius norm")
            self.__comp_CFN()
            self.get_ordered_cfn() ### Let's always do it
            
        del self.__invC ### TODO: maybe sometimes one wants to save it?
        del self.__Pi   ### 

            
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
            
    def __add_pseudocounts(self):
        """Adds pseudocounts to the observed mutation frequencies"""
        self.__Pi = (1.-self.pseudocount_weight)*self.__Pi +\
                  self.pseudocount_weight/self.q*np.ones((self.N,self.q))
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
        self.__invC=tensorinv(C) # THIS IS MINUS THE COUPLINGS J_ij
        self.gauge='Gas-lattice' ### we are in the Gas-lattice gauge

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
          
    def __comp_DI(self):
        """Computes Direct Information"""
        ### TODO: implement other DCAmethods
        ### TODO: there should be a smarter way of storing and accessing these symmetric matrices using only half the space
        self.direct_information=np.zeros((self.N,self.N))
        ### TODO: this loops are not very python friendly. There might be a fastest way...
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                # direct information from mean-field
                self.direct_information[i,j] = self.__bp_link(i,j)
        self.direct_information+=self.direct_information.T

    def __to_ising_gauge(self):
        """converts the couplings to the Ising gauge"""
        ### TODO: how do I go back to the lattice-gas Gauge?
        # Remember: invC=-J_ij
        # invC is symmetric in the mfDCA!
        J_avg=np.average(self.__invC,axis=3)
        ### Possa dio aver pieta' di noi
        self.__invC+=np.average(J_avg,axis=1)[:,np.newaxis,:,np.newaxis]\
                      -J_avg[:,:,:,np.newaxis]\
                      -np.transpose(J_avg,axes=(0,2,1))[:,np.newaxis,:,:]
        self.gauge='Ising'
                      
        
    def __comp_CFN(self,no_gaps=False):
        """Computes Frobenius norm"""
        ### TODO: add option to exclude gaps
        ### TODO: CFN can be negative? Check this!
        if self.gauge!='Ising':
            self.__to_ising_gauge()
        self.CFN=np.sqrt(np.sum(self.__invC**2,axis=(1,3))) # NB: invC^2 so the sign doesn't matter
        ### TODO: is CFN symmetric??
        F_sum=np.sum(self.CFN,axis=0)
        self.CFN-=F_sum[:,np.newaxis]*F_sum[np.newaxis,:]/np.sum(F_sum)
        
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
        ###      apparently there is a mathematical reason for it
        DI = np.trace(\
                      np.dot(Pdir.T , np.log((Pdir+tiny)/(Pfac+tiny)) ) \
        )
        return DI

    def get_ordered_di(self,k_pairs=None,offdiag=4,return_di=False):
        """Sort the pairs by their direct information"""
        faraway=np.triu_indices(self.N,k=offdiag)
        self.di_order=(np.array(faraway).T[np.argsort(self.direct_information[faraway])])[::-1]
        if return_di:
            if k_pairs==None:
                k_pairs=self.N*2
            return self.direct_information[[self.di_order[:k_pairs,0],self.di_order[:k_pairs,1]]] ### TODO: this is not creating a copy. Be careful

    def get_ordered_cfn(self,k_pairs=None,offdiag=4,return_score=False):
        """Sort the pairs by their Frobenius norm score"""
        faraway=np.triu_indices(self.N,k=offdiag)
        self.cfn_order=(np.array(faraway).T[np.argsort(self.CFN[faraway])])[::-1]
        if return_score:
            if k_pairs==None:
                k_pairs=self.N*2
            return self.CFN[[self.cfn_order[:k_pairs,0],self.cfn_order[:k_pairs,1]]] ### TODO: this is not creating a copy. Be careful
        
    def print_results(self,filename):
        """Prints DI and MI (compatible with the output of the matlab code)"""
        fh=open(filename,'w')
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                fh.write('%d %d'%(i+1, j+1)) # matlab indexing    
                fh.write(' %g'%self.mutual_information[i,j])
                fh.write(' %g'%self.direct_information[i,j])
                fh.write('\n')
        fh.close()

    def print_DI(self,filename):
        """Prints DI (compatible with SPECTRUS-evo input)"""
        fh=open(filename,'w')
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                fh.write('%d %d'%(i, j))                
                fh.write(' %g'%self.direct_information[i,j])
                fh.write('\n')
        fh.close()

    def print_contacts(self,filename,iseq,n_pairs=None,score='DI',):
        """Prints pairs with highest coupling score (compatible with AWSEM-ER input)"""
        
        if n_pairs==None:
            n_pairs=self.N*2
        ### I'm not sure this is the most elegant and correct thing to do to check/select the score option
        if score=='DI':
            ix=self.di_order[:n_pairs,0]
            iy=self.di_order[:n_pairs,1]
        elif score=='CFN':
            ix=self.cfn_order[:n_pairs,0]
            iy=self.cfn_order[:n_pairs,1]
        else:
            raise ValueError("Unrecognize score option: '%s'"%score)
        ###
        
        assert iseq>=0 and iseq<self.M,'ERROR: invalid sequence ID'
    
        fh=open(filename,'w')
        fh.write('i   j   i_id  j_id\n')
        for i,j in zip(ix,iy):
            i0=self.alignment.align2orig[iseq][self.alignment.strip2align[i]]
            j0=self.alignment.align2orig[iseq][self.alignment.strip2align[j]]
            res_i=self.alignment.stripped_seqs[iseq][i]
            res_j=self.alignment.stripped_seqs[iseq][j]
            fh.write("%d %d %d_%s %d_%s\n"%\
                     (i0,j0,i0,res_i,j0,res_j))
        fh.close()

        
def plot_contacts(dca_obj,n_pairs=None,lower_half=False,iseq=None,colormap=plt.cm.CMRmap_r,binary=False,offset=0,score='DI'):
    """Prints the contact maps derived from a DCA object.

    if iseq>0 will remap the indexes to the original aminoacids of the sequence;

    if ofset>0 will shift the index of the first aminoacid (use it to compare dca on different part of a sequence);

    if you want to compare the contacts from two dca objects, just use

        plot_contacts(dca_obj1)
        plot_contacts(dca_obj2,lower_half=True)

    score options:
    'DI' (default) -> direct information as defined by Morcos et al., PNAS 2011
    'CFN' -> Frobenius norm as defined by Ekerberg et al., PRE 2013

"""
    ### TODO: how can we change this to plot and compare two contact maps?
    ###       is it better to do it inside the function or outside?
    ### TODO: be careful with remapping indeces!
    ###       right now if there are "-" in the sequence considered
    ###       the indeces are actually counted twice!!
    if n_pairs==None:
        n_pairs=dca_obj.N*2
    ### I'm not sure this is the most elegant and correct thing to do to check/select the score option
    if score=='DI':
        ix=dca_obj.di_order[:n_pairs,0]
        iy=dca_obj.di_order[:n_pairs,1]
    elif score=='CFN':
        ix=dca_obj.cfn_order[:n_pairs,0]
        iy=dca_obj.cfn_order[:n_pairs,1]
    else:
        raise ValueError("Unrecognize score option: '%s'"%score)
    ###
    old_ix=ix
    old_iy=iy
    if iseq!=None:
        assert iseq>=0 and iseq<dca_obj.M,'ERROR: invalid sequence ID'
        ix=dca_obj.alignment.align2orig[iseq][dca_obj.alignment.strip2align[ix]]#+offset
        iy=dca_obj.alignment.align2orig[iseq][dca_obj.alignment.strip2align[iy]]#+offset
    matr=np.zeros((dca_obj.N,dca_obj.N)) ### TODO: here should not be N but the length of the original sequence
    if lower_half:
        if binary:
            matr[[ix,iy]]=1
        else:
            if score=='DI':
                matr[[ix,iy]]=dca_obj.direct_information[[old_ix,old_iy]]
            if score=='CFN':
                matr[[ix,iy]]=dca_obj.CFN[[old_ix,old_iy]]

        iny, inx = np.indices(matr.shape) 
        my_mask=inx<=iny
    else:
        if binary:
            matr[[iy,ix]]=1
        else:
            if score=='DI':
                matr[[iy,ix]]=dca_obj.direct_information[[old_iy,old_ix]]
            if score=='CFN':
                matr[[iy,ix]]=dca_obj.CFN[[old_iy,old_ix]]
        iny, inx = np.indices(matr.shape) 
        my_mask=inx>=iny
        #plt.scatter(ix,iy,marker='s',s=3,color=colore)
    plt.imshow(np.ma.array(matr,mask=my_mask),cmap=colormap,origin='lower',extent=[offset,dca_obj.N+offset,offset,dca_obj.N+offset]) ### TODO: here should not be N but the length of the original sequence
    plt.plot(range(dca_obj.N),color='black') ### TODO: here should not be N but the length of the original sequence
    return matr


def compute_dca(inputfile,pseudocount_weight=0.5,theta=0.1,compute_MI=False,compute_CFN=False):
    """Perform mfDCA starting from a FASTA input file. Returns a DCA object"""
    ### TODO: add "filter" argument to filter sequences with too many gaps. add "method" arguments to use different DCA implementations
    alignment=sf.read_alignment(inputfile,check_aminoacid=True) ### TODO: add filter_limit here
    print("=== DCA analysis ===\n Number of sequences = %d\n Alignment lenght= %d\n"%(alignment.M,alignment.N))
    dca_obj=DCA(alignment,pseudocount_weight=pseudocount_weight,theta=theta,get_MI=compute_MI,get_DI=True,get_CFN=compute_CFN)
    print(" Effective number of sequences = %d\n"%dca_obj.Meff)
    #dca_obj.get_ordered_di()
    print("=== DCA completed ===")
    return dca_obj
