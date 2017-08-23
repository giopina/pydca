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
	def __init__(self, fasta_filename, structure_filename, alignment_filename, score_filename, muscle_path, output_filename="correlation.txt"):
		"""
		Input formats:
		- FASTA file: Must contain only the sequence of the chosen chain, in FASTA format (ungapped, as downloaded from the PDB)
		- Structure file: .pdb file containing ONLY the chain of choice, as downloaded from the PDB
		- Alignment file: Must contain only the sequence of choice, in FASTA format (gapped, as extracted from its Pfam family)
		- Score file: Must have 3 columns: int, int, float. The first two are the indices of the score matrix, the third are the score values.
		              All off-diagonal entries must be reported. Order can be arbitrary (e.g., sorted by score...)
		Note: you have to specify the path for using the MUSCLE alignment program
		Output format:
		- Output: Text file with two columns. The first reports structural distances,, the second reports the corresponding DCA scores.
		"""
		self.alignment_filename = alignment_filename    #'5p21_chosen.fasta.txt'
		self.score_filename = score_filename    #'formatted_score.txt'
		self.structure_filename = structure_filename    #'5p21.pdb'
		self.fasta_filename = fasta_filename    #'5p21.fasta.txt'
		self.output_filename = output_filename     #'check.out'
		self.muscle_path = muscle_path    #'/Users/edoardosarti/programs/muscle/muscle3.8.31_i86darwin64'

		self.io = PDBIO()
		self.ppb = PPBuilder()

		DI_score_matrix, contact_matrix = self.calculate_correlation()


	def muscle_alignment(self, sequences, workdir, cleanup=False):
		"""
		Input:
		- list(strings)::sequences list of sequences
		- string::workdir path of the working directory where a temp folder will be created in order to run the program
		Output:
		- list(strings)::seqaln list of aligned sequences, in the same order of the input
		"""

		tmpdir = workdir + '/tmp_seqaln/'
		tmpin_filename = tmpdir + 'seqs_in.fa'
		tmpout_filename = tmpdir + 'seqs_out.fa'

		# Create the temp folder and write the input file
		if os.path.exists(tmpdir):
			shutil.rmtree(tmpdir)
		os.mkdir(tmpdir)
		tmpin_file = open(tmpin_filename, 'w')
		for ns in range(len(sequences)):
			tmpin_file.write('>Sequence {0}\n'.format(ns))
			tmpin_file.write('{0}\n'.format(sequences[ns].upper()))
		tmpin_file.close()

		tries = 0
		while (not os.path.exists(tmpout_filename)):
			if tries >= 30:
				raise NameError("There is something wrong with this: {0}".format(tmpin_filename))
			fnull = open(os.devnull, 'w')
			p = subprocess.Popen([self.muscle_path, '-in', tmpin_filename, '-out', tmpout_filename], stdout=fnull, stderr=fnull)
			p.wait()
			fnull.close()
			tries += 1

		tmpout_file = open(tmpout_filename, 'r')
		text = tmpout_file.read().split('\n')
		tmpout_file.close()

		if cleanup:
			os.rmdir(tmpdir)

		seqaln = []
		for line in text:
			if not line:
				continue
			if line[0] == '>':
				seqaln.append('')
			else:
				seqaln[-1] += line.strip()

		return seqaln

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


	def parse_sequence(self, sequence_filename, tag=None, from_MSA=False):
		print(sequence_filename)
		sequence_file = open(sequence_filename, 'r')
		text = sequence_file.read().split('\n')
		sequence = ''
		is_tagged = False
		for line in text:
			if not line:
				continue
			fields = line.split()
			if fields[0][0] == '>':
				if tag != None and tag in fields[0]:
					is_tagged = True
				else:
					# If the name does not include the tag and the previous name did, the search is over
					if is_tagged:
						break
					continue
			if tag == None or (tag != None and is_tagged):
				sequence += fields[0].strip()

		if from_MSA:
			sequence = sequence.replace('.','').replace('-','').replace('*','')

		#print(sequence)
		return sequence


	def calculate_correlation(self):
		# Get sequence from MSA
		# The MSA only contains part of the whole sequence
		sequence_from_msa = self.parse_sequence(self.alignment_filename)
		sequence_from_msa_for_aln = sequence_from_msa.replace('.','').replace('-','')

		# Get sequence from FASTA file
		# This is the only reference sequence
		sequence_from_fasta = self.parse_sequence(self.fasta_filename)

		aligned_sequences = self.muscle_alignment([sequence_from_fasta, sequence_from_msa_for_aln], os.getcwd())
		#print(aligned_sequences)

		sequence_from_msa_for_pairing = [x for x in sequence_from_msa.replace('.','')][::-1]
		nMSA = 0    # Index of MSA positions: capital letters + dash-gaps
		nrMSA = 0    # Index of MSA residues: capital letters only
		nPDB = 0    # Index of PDB residues (sequential one, not the true resid)
		MSAgaps = []    # List of MSA indices where dash-gaps are found
		from_MSA_to_PDB = {}    # From MSA to PDB is not a bijective function, since MSA has '-' positions, which are of course overlooked in the PDB
		from_rMSA_to_PDB = {}    # This is a bijective function
		from_PDB_to_rMSA = {}    # This is the inverse function
		for nc in range(len(aligned_sequences[1])):
		#print(aligned_sequences[0][nc], aligned_sequences[1][nc])
			if aligned_sequences[0][nc] != '-':
				nPDB += 1
			if aligned_sequences[1][nc] != '-':    # The dashes contained in the muscle alignment are not the dashes reported on the MSA ;)
				found = False
				while(not found):
					cpop = sequence_from_msa_for_pairing.pop()    # Contains [A-Za-z\-] characters
					#print(cpop)
					if cpop == aligned_sequences[1][nc] and re.match('[A-Z]',aligned_sequences[1][nc]):
						from_MSA_to_PDB[nMSA] = nPDB
                                                from_rMSA_to_PDB[nrMSA] = nPDB
                                                from_PDB_to_rMSA[nPDB] = nrMSA
						nMSA += 1
						nrMSA += 1
						found = True
					elif cpop == '-':
						from_MSA_to_PDB[nMSA] = -1
						MSAgaps.append(nMSA)
						nMSA += 1
					elif cpop.upper() == aligned_sequences[1][nc]:
						found = True

		# Get sequence from PDB file, via Biopython
		# The structural sequence can have holes and is not suitable for being aligned
		parser = PDBParser(PERMISSIVE=True)
		structure = parser.get_structure('self', self.structure_filename)
		for c in structure[0]:
			chain_id = c.get_id()
			break

		sequence_from_pdb = ''    # In the PDB file there must be only one model with only one chain: the one of interest. But this chain could contain holes
		from_struct_to_PDB = {}
		from_PDB_to_struct = {}
		resids_in_MSA = []
		nres = 0
		for ppept in self.ppb.build_peptides(structure[0][chain_id]):
			sequence_from_pdb += ppept.get_sequence()
			for res in ppept:
				resid = res.get_id()[1]
				from_struct_to_PDB[nres] = resid
				from_PDB_to_struct[resid] = nres
				if resid in list(from_MSA_to_PDB.values()):
					resids_in_MSA.append(resid)
				nres += 1
	
		structure_D = self.min_distance_method(structure[0][chain_id])

		DI_list = []
		score_file = open(self.score_filename, 'r')
		text = score_file.read().split('\n')
		for line in text:
			if not line:
				continue
			fields = line.split()
			DI_list.append((int(fields[0]), int(fields[1]), float(fields[2])))

		last_index = sorted(DI_list, key = lambda x : x[1])[-1][1]
		print("Input score matrix (from MSA) has dimension:", last_index)
		print("Number of MSA residues in the chosen sequence:", len(resids_in_MSA))
		DI = np.zeros((last_index, last_index))
		for entry in sorted(DI_list, key = lambda x : (x[0],x[1])):
			DI[entry[0]-1][entry[1]-1] = entry[2]
			DI[entry[1]-1][entry[0]-1] = entry[2]

		rDI = np.array(DI)
		for g in MSAgaps[::-1]:
			rDI = np.delete(rDI, g, 0)
			rDI = np.delete(rDI, g, 1)
		#print(rDI.shape)

		for r in range(len(rDI[0])-1,-1,-1):
			if from_rMSA_to_PDB[r] not in resids_in_MSA:
				rDI = np.delete(rDI, r, 0)
				rDI = np.delete(rDI, r, 1)
		print("Number of sequence-mapped MSA residues mapped in the structure:", rDI.shape[0])

		S = np.zeros(rDI.shape)
		nr1 = 0
		from_PDB_to_S = {}
		from_S_to_PDB = {}
		for r1 in sorted(resids_in_MSA):
			nr2 = 0
			for r2 in sorted(resids_in_MSA):
				S[nr1][nr2] = structure_D[from_PDB_to_struct[r1]][from_PDB_to_struct[r2]]
				nr2 += 1
			from_PDB_to_S[r1] = nr1
			from_S_to_PDB[nr1] = r1
			nr1 += 1

		#print(S)
		#print(rDI)

		print("Correlation file:", self.output_filename)
		check_file = open(self.output_filename, 'w')
		for nr1 in range(len(S[0])):
			for nr2 in range(len(S[0])):
				check_file.write("{0:.6f}\t{1:.6f}\n".format(S[nr1][nr2], rDI[nr1][nr2]))
		check_file.close()

		return S, rDI
