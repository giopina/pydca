import mdtraj
import numpy as np
#import dca

class pdb_contacts:

    def __init__(self,pdb_name):
        # read the structure from a pdb file (one chain, no hetatm, no water, just the good old protein, please
        self.struct=mdtraj.load(pdb_name)
        # compute the contacts. This computes all the minimum atom distances between residues. The cutoff will be applied later
        self.dd,self.rp=mdtraj.compute_contacts(self.struct)
        # ah si, this is the sequence of the pdb chain (the residues which are resolved, may be less than the FASTA from DB website)
        self.seq=''.join([r.code for r in self.struct.topology.residues])
        # not using this but maybe someday can be useful
        #cc=mdtraj.geometry.squareform(dd,rp)

    def get_pair_idx(self,cutoff=0.8):
        return self.rp[self.dd[0]<cutoff,0],self.rp[self.dd[0]<cutoff,1]


    
def compareDCAandPDB(dca_obj,pdb_obj,n_pairs=None,iseq=0,score='DI',cutoff=0.5,pdb_seq_subset=None):
    """Function to compare contacts. This uses a scatterplot and also compares the two contact maps and find the intersection.                                                                                
    """
    ix,iy,old_ix,old_iy=dca_obj.get_pair_idx(n_pairs=n_pairs,iseq=iseq,score=score)
    idx1=np.array((ix,iy)).T
    ix,iy=pdb_obj.get_pair_idx(cutoff=cutoff)
    if pdb_seq_subset != None:
        # defining the remapping array
        remapp=np.ones(len(pdb_obj.seq),dtype=np.int32)*-1
        remapp[pdb_seq_subset]=np.arange(pdb_seq_subset.shape[0],dtype=np.int32)
        # remapping indexes on the seq subset of interest
        ix=remapp[ix]
        iy=remapp[iy]
        # removing negative indeces
        iy=iy[ix>=0]
        ix=ix[ix>=0]
        ix=ix[iy>=0]
        iy=iy[iy>=0]
    else:
        ix_new=ix
        iy_new=iy
    idx2=np.array((ix,iy)).T
    idx_both=np.array([x for x in set(tuple(x) for x in idx1) & set(tuple(x) for x in idx2)])
    n_common=idx_both.shape[0]
    TPR=idx_both.shape[0]/len(idx1)
    return TPR, n_common
