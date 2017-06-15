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

