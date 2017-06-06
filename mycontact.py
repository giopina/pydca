import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import shutil
import re
from Bio.PDB import *
import warnings
from Bio.PDB.PDBParser import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)


class Compare:
	def __init__(self, fasta_filename, output_filename="correlation.txt"):
		"""
		Input formats:
		- FASTA file: Must contain only the sequence of the chosen chain, in FASTA format (ungapped, as downloaded from the PDB)
		- Structure file: .pdb file containing ONLY the chain of choice, as downloaded from the PDB
		Output format:
		- Output: Text file with two columns. The first reports structural distances,, the second reports the corresponding DCA scores.
		"""
		self.structure_filename = structure_filename    #'5p21.pdb'
		self.fasta_filename = fasta_filename    #'5p21.fasta.txt'
		self.output_filename = output_filename     #'check.out'

		self.io = PDBIO() # was is das?
		self.ppb = PPBuilder()

		DI_score_matrix, contact_matrix = self.calculate_correlation()

	def calculate_correlation(self):

		# Get sequence from FASTA file
		# This is the only reference sequence
		sequence_from_fasta = self.parse_sequence(self.fasta_filename) ###

		# Get sequence from PDB file, via Biopython
		# The structural sequence can have holes and is not suitable for being aligned
		parser = PDBParser(PERMISSIVE=True)
		structure = parser.get_structure('self', self.structure_filename)
		for c in structure[0]:
			chain_id = c.get_id()
			break

		#sequence_from_pdb = ''    # In the PDB file there must be only one model with only one chain: the one of interest. But this chain could contain holes
		#from_struct_to_PDB = {}
		#from_PDB_to_struct = {}
		#resids_in_MSA = []
		#nres = 0
		#for ppept in self.ppb.build_peptides(structure[0][chain_id]):
		#	sequence_from_pdb += ppept.get_sequence()
		#	for res in ppept:
		#		resid = res.get_id()[1]
		#		from_struct_to_PDB[nres] = resid
		#		from_PDB_to_struct[resid] = nres
		#		if resid in list(from_MSA_to_PDB.values()):
		#			resids_in_MSA.append(resid)
		#		nres += 1
	
		structure_D = self.min_distance_method(structure[0][chain_id])


		return S, rDI

        
        def CA_distance_method(self, bp_chain):
		"""
		Calculates distances between CAs.
		Input:
		- biopython_chain::bp_chain chain object from biopython
		Output:
		- numpy.array((N,N))::D distance matrix between the N residues in the chain
		"""

		rescoords = []
		for ppept in self.ppb.build_peptides(bp_chain): 
			for residue in ppept:
				rescoords.append(residue['CA'].get_vector().get_array())
		N = len(rescoords)
		D = np.zeros((N,N))
		for nr1 in range(len(rescoords)):
			for nr2 in range(nr1+1, len(rescoords)):
				d = np.linalg.norm(rescoords[nr1] - rescoords[nr2])
				D[nr1][nr2] = d
				D[nr2][nr1] = d
		return D

	def min_distance_method(self, bp_chain):
		"""
		Calculates the minimum distance between side chains.
		Input:
		- biopython_chain::bp_chain chain object from biopython
		Output:
		- numpy.array((N,N))::D distance matrix between the N residues in the chain
		"""

		backbone = ['CA', 'N', 'O', 'H']
		rescoords = []
		for ppept in self.ppb.build_peptides(bp_chain):
			for residue in ppept:
				rescoords.append(residue['CA'].get_vector().get_array())
		N = len(rescoords)
		D = np.zeros((N,N))
		nr1 = 0
		for ppept in self.ppb.build_peptides(bp_chain):
			for residue in ppept:
				nr2 = 0
				for ppept2 in self.ppb.build_peptides(bp_chain):
					for residue2 in ppept2:
						if residue == residue2: #or np.linalg.norm(residue['CA'].get_vector().get_array() - residue2['CA'].get_vector().get_array()) > 15:
							continue
						dmin = 10000
						for atom in residue:
							if atom.get_id() in backbone:
								continue
							for atom2 in residue2:
								if atom2.get_id() in backbone:
									continue
								d = np.linalg.norm(atom.get_vector().get_array() - atom2.get_vector().get_array())
								if d < dmin:
									dmin = d
						D[nr1][nr2] = dmin
						D[nr2][nr1] = dmin
						nr2 += 1
				nr1 += 1
		return D		
