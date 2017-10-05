import pydca.dca as _dca
import pydca.mycontact as _mc
import matplotlib.pyplot as _plt
import numpy as _np
        
def plot_contacts(dca_obj,n_pairs=None,lower_half=False,iseq=None,colormap=_plt.cm.CMRmap_r,binary=True,offset=0,score='DI'):
    """Prints the contact maps derived from a DCA object.

    if iseq>0 will remap the indexes to the original aminoacids of the sequence;

    if offset>0 will shift the index of the first aminoacid (use it to compare dca on different part of a sequence);

    if you want to compare the contacts from two dca objects, just use

        plot_contacts(dca_obj1)
        plot_contacts(dca_obj2,lower_half=True)

    score options:
    'DI' (default) -> direct information as defined by Morcos et al., PNAS 2011
    'CFN' -> Frobenius norm as defined by Ekerberg et al., PRE 2013

    lower_half=True prints the contact map in the bottom-right triangle of the plot,
    default prints it on top-left side
"""
    ### TODO: how can we change this to plot and compare two contact maps?
    ###       is it better to do it inside the function or outside?
    ix,iy,old_ix,old_iy=dca_obj.get_pair_idx(n_pairs=n_pairs,iseq=iseq,score=score)
    matr=_np.zeros((dca_obj.alignment.N_orig[iseq],dca_obj.alignment.N_orig[iseq])) # matrix of zeros

    if lower_half:
        if binary:
            matr[[ix,iy]]=1 # black and white plot
        else:
            # color based on the score
            if score=='DI':
                matr[[ix,iy]]=dca_obj.direct_information[[old_ix,old_iy]]
            if score=='CFN':
                matr[[ix,iy]]=dca_obj.CFN[[old_ix,old_iy]]
        iny, inx = _np.indices(matr.shape) 
        my_mask=inx<=iny # This will fill only the desired half of the canvas
    else:
        if binary:
            matr[[iy,ix]]=1 # black and white plot
        else:
            # color based on the score
            if score=='DI':
                matr[[iy,ix]]=dca_obj.direct_information[[old_iy,old_ix]]
            if score=='CFN':
                matr[[iy,ix]]=dca_obj.CFN[[old_iy,old_ix]]
        iny, inx = _np.indices(matr.shape) 
        my_mask=inx>=iny # This will fill only the desired half of the canvas

    # Now let's plot the contact map...
    _plt.imshow(_np.ma.array(matr,mask=my_mask),cmap=colormap,origin='lower',\
               extent=[offset,dca_obj.alignment.N_orig[iseq]+offset,\
                       offset,dca_obj.alignment.N_orig[iseq]+offset])
    # ...and we draw a line on the diagonal, just for fun
    _plt.plot(range(dca_obj.alignment.N_orig[iseq]),color='black')
    return matr # return matrix of contacts if one wants to replot it differently

def scatter_contacts(dca_obj1,dca_obj2,n_pairs=(None,None),iseq=(0,0),score='DI'):
    """Another function to plot dca/pdb contacts. This uses a scatterplot and also compare the two contact maps and find the intersection
    """
    ### TODO: this can be modified to be used also with a contact map from a PDB structure.
    ix,iy,old_ix,old_iy=dca_obj1.get_pair_idx(n_pairs=n_pairs[0],iseq=iseq[0],score=score)
    idx1=_np.array((ix,iy)).T
    ix,iy,old_ix,old_iy=dca_obj2.get_pair_idx(n_pairs=n_pairs[1],iseq=iseq[1],score=score)
    idx2=_np.array((ix,iy)).T
    idx_both=_np.array([x for x in set(tuple(x) for x in idx1) & set(tuple(x) for x in idx2)])
    print('Number common contacts = %d'%(idx_both.shape[0]))
    print('Fraction of common contacts = %.2f'%(idx_both.shape[0]/(len(idx1)+len(idx2))*2))
    _plt.figure(figsize=(10,10))
    _plt.scatter(idx1[:,0],idx1[:,1],alpha=0.99,s=10,c='cyan',marker='s',label='DCA 1')
    _plt.scatter(idx2[:,1],idx2[:,0],alpha=0.99,s=10,c='green',marker='s',label='DCA 2')
    _plt.scatter(idx_both[:,0],idx_both[:,1],alpha=0.99,s=18,c='red',label='common contacts')
    _plt.scatter(idx_both[:,1],idx_both[:,0],alpha=0.99,s=18,c='red')
    _plt.legend()

def plotDCAandPDB(dca_obj,pdb_obj,n_pairs=None,iseq=0,score='DI',cutoff=0.5,pdb_seq_subset=None):
    """Function to plot dca/pdb contacts. This uses a scatterplot and also compares the two contact maps and find the intersection.                                                                               
    """
    ix,iy,old_ix,old_iy=dca_obj.get_pair_idx(n_pairs=n_pairs,iseq=iseq,score=score)
    idx1=_np.array((ix,iy)).T
    ix,iy=pdb_obj.get_pair_idx(cutoff=cutoff)
    if pdb_seq_subset != None:
        # defining the remapping array
        remapp=_np.ones(len(pdb_obj.seq),dtype=_np.int32)*-1
        remapp[pdb_seq_subset]=_np.arange(pdb_seq_subset.shape[0],dtype=_np.int32)
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
    idx2=_np.array((ix,iy)).T
    idx_both=_np.array([x for x in set(tuple(x) for x in idx1) & set(tuple(x) for x in idx2)])
    ### TODO: up to here this is the same as mycontacts.py->compareDCAandPDB, so maybe they can be merged together somehow...
    print( 'Number common contacts = %d' % (idx_both.shape[0]) )
    print( 'Common contacts / DCA pairs = %.2f' % \
           (idx_both.shape[0]/len(idx1)) )
    _plt.figure(figsize=(10,10))
    _plt.scatter(idx1[:,0],idx1[:,1],alpha=0.99,s=10,c='cyan',marker='s',label='DCA')
    _plt.scatter(idx2[:,1],idx2[:,0],alpha=0.99,s=10,c='green',marker='s',label='PDB')
    _plt.scatter(idx_both[:,0],idx_both[:,1],alpha=0.99,s=18,c='red',label='common contacts')
    _plt.scatter(idx_both[:,1],idx_both[:,0],alpha=0.99,s=18,c='red')
    _plt.legend()
