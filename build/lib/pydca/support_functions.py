import re
import datetime
import numpy as _np

def write_log(name, text):
        """ Quietly writes a log """
        log_filename = name + '.log'
        log_file = open(log_filename, 'w')
        log_file.write(text)
        log_file.close()
        

def header(name):
	""" Creates a tidy header for logs, with date and time"""
	return "[" + name + " " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "] "


def raise_error(name, text):
	""" Version of write_log for reporting errors and stopping"""
	indent = " "*len(header(name))
	lines = text.split('\n')
	logmsg = header(name) + lines[0] + '\n'
	for line in lines[1:]:
		logmsg += indent + line + '\n'
	write_log(name, logmsg)
	raise NameError(logmsg)


def print_log(name, text, quiet=False):
	""" Handy function for writing logs on screen and in a log file.
	You have to give a name to the function you're logging from"""
	indent = " "*len(header(name))
	lines = text.split('\n')
	logmsg = header(name) + lines[0]
	for line in lines[1:]:
		logmsg += '\n' + indent + line
	if not quiet:
		print(logmsg)
	return logmsg + '\n'


def FASTA_parser(fasta_filename, to_upper=False, check_aminoacid=False, check_nucleicacid=False):
	""" Basic FASTA parser. 
	Input: path of a fasta-formatted file
	Output: list of dictionaries with keys 'title' and 'sequence'
	Optional flags:
	    to_upper: writes sequence in uppercase
	    check_aminoacid: checks if the sequence is composed by FASTA-compliant aminoacid characters
	    check_nucleicacid: checks if the sequence is composed by FASTA-compliant nucleic acid characters
	"""
	this_name = "FASTA_parser"

	# check_aminoacid and check_nucleicacid are mutually exclusive
	if check_aminoacid and check_nucleicacid:
		raise_error(this_name, "ERROR: only one out of 'check_aminoacid' and 'check_nucleicacid' options can be active")

	fasta_list = []   # The list of dictionaries (output)

	# Matches entire non-empty lines composed only by alphabetical characters, hyphens,
	# asterisks, backslashes and dots. Does not match FASTA titles!
	fasta_pattern = re.compile('^[A-Za-z\-\*\.\\\\]+$')

	fasta_file = open(fasta_filename, 'r')
	text = fasta_file.read().split('\n')
	fasta_file.close()

	there_is_title_line = False
	nl = 0
	no_titles = 0
	for t_line in text:
		nl += 1
		line = t_line.strip()
		# Standard check for empty lines
		if not line:
			there_is_title_line = False
			continue
		# If it is not a FASTA line or a FASTA title, stop
		if not (line[0] == '>' or re.match(fasta_pattern, line)):
			raise_error(this_name, "ERROR: Line {0} does not respect the FASTA format".format(nl))
		# If it is a FASTA title, expect only one title, and put the title in the 'title' key of a new entry dict in the main list
		# If there is no title (just a ">"), make up one
		if line[0] == '>':
			if there_is_title_line == False:
				if len(line) > 1:
					title = line[1:]
				else:
					no_titles += 1
					title = 'NO_TITLE_' + str(no_titles)
				fasta_list.append({'title' : title, 'sequence' : ''})
				there_is_title_line = True
			else:
				raise_error(this_name, "ERROR: sequence title with no amino acid sequence at line {0}".format(nl))
		# If it is a FASTA line, join that line to the 'sequence' key of the latest entry dict in the main list
		elif re.match(fasta_pattern, line):
			if not fasta_list or (fasta_list[-1]['sequence'] == '' and not there_is_title_line):
				raise_error(this_name, "ERROR: there is a string of text not preceded by a title at line {0}".format(nl))
			if to_upper:
				fasta_list[-1]['sequence'] += line.strip().replace('\\', '').upper()
			else:
				fasta_list[-1]['sequence'] += line.strip().replace('\\', '')
			there_is_title_line = False

	# Optional format checks
	allowed_characters = None
	if check_aminoacid:
		allowed_characters = 'ABCDEFGHIKLMNPQRSTUVWYZXabcdefghiklmnpqrstuvwyzx*.-'
	elif check_nucleicacid:
		allowed_characters = 'ACGTURYKMSWBDHVNacgturykmswbdhvn-'
	if allowed_characters:
		for seq in fasta_list:
			for lett in seq['sequence']:
				if lett not in allowed_characters:
					print_log(this_name, "WARNING: Sequence with title {0} contains unknown character {1}".format(seq['title'], lett))
	return fasta_list


class Alignment:
        
        def __init__(self,fasta_list):
                """This creates an Alignment object from a list of sequences. It strip them from lowercase letters and '.' or '*' characters.
                -----------------------
                A bit of nomenclature:

                aligned sequence -> sequence in the alignment fasta file (with "-","." and lowercase letters)
                original sequence -> the original sequence of the protein (without "-" and "." gaps, but with lowercase letters)
                stripped sequence -> sequence stripped from "." and lowercase letters. Ready for the DCA
                -----------------------
                Some variables:
                M        -> number of sequences
                N        -> length of stripped sequences
                N_align  -> length of sequences in the alignment
                N_orig   -> list of lengths of original sequences
                """
                self.letter2numer={\
                # full AA alphabet
                # TODO: check all the possible letters not in the 20 aa... Maybe this is just stupid...
                                   '-':0 ,\
                                   'A':1 ,\
                                   'C':2 ,\
                                   'D':3 ,\
                                   'E':4 ,\
                                   'F':5 ,\
                                   'G':6 ,\
                                   'H':7 ,\
                                   'I':8 ,\
                                   'K':9 ,\
                                   'L':10,\
                                   'M':11,\
                                   'N':12,\
                                   'P':13,\
                                   'Q':14,\
                                   'R':15,\
                                   'S':16,\
                                   'T':17,\
                                   'V':18,\
                                   'W':19,\
                                   'Y':20,\
                                   'X':0,\
                                   'Z':0,\
                                   'B':0,\
                }
                self.numer2letter=list(self.letter2numer.keys())[:-1] # NB: 0 will always be backmapped to "-"
                
                self.M=len(fasta_list) # this is the number of sequences
                
                self.sequences=[seq['sequence'] for seq in fasta_list] # sequences
                self.names=[seq['title'] for seq in fasta_list] # names of sequences can be useful in the prostprocessing

                if len(set([len(s) for s in self.sequences]))>1:
                        raise_error(this_name,"ERROR: aligned sequences have different lengths!")
                self.N_align=len(self.sequences[0]) # this is the lenght of each aligned sequences

                self.orig2align=[_np.where([i!='-' and i!='.' for i in stringa])[0] for stringa in self.sequences] # this takes x3 times than stripping the seq.
                ### the following def tranform gaps in the index of the previous non-gap (not a very good choice)
                #self.align2orig=[_np.cumsum([i!='-' and i!='.' for i in stringa])-1 for stringa in self.sequences] # ok, this is likely to be unefficient but who cares...
                ### the following instead transform gaps into -1s (hopefully)
                self.align2orig=[_np.ones(self.N_align,dtype=_np.int32)*-1 for i in range(self.M)]
                self.N_orig=[]
                #for iseq,seq in enumerate(self.sequences):
                for iseq in range(self.M):
                        #n_orig=_np.sum([i!='-' and i!='.' for i in seq])
                        n_orig=len(self.orig2align[iseq])
                        self.N_orig.append(n_orig)
                        #self.align2orig[iseq][_np.where([i!='-' and i!='.' for i in seq])[0]]=_np.arange(n_orig)
                        self.align2orig[iseq][self.orig2align[iseq]]=_np.arange(n_orig)
                
                self.__strip()
                
        def __strip(self):
                """Fuck this. I'm going to be a stripper"""
                this_name='stripper'
                self.stripped_seqs=[re.sub("[a-z]|\.|\*",'',s) for s in self.sequences] # this should remove lowercase letters and "."  and "*"
                ### TODO: do we really need to do search for . and * for all sequences??
                ### Maybe it's better to do it, just in case there's some error in the fasta file.
                ### It doesn't take too much time anyway
                if len(set([len(s) for s in self.stripped_seqs]))>1:
                        raise_error(this_name,"ERROR: stripped sequences have different lengths!")
                self.N=len(self.stripped_seqs[0]) # this is the lenght of each stripped sequences

                ### these two lists of lists are redundant. Since the gaps are always in the same positions only one list for all the sequences is enough...
                #self.strip2align=[_np.where([i!='.' and i!='*' and not i.islower() for i in stringa])[0] for stringa in self.sequences] # this takes x3 times than stripping the seq.
                #self.align2strip=[_np.cumsum([i!='.' and i!='*' and not i.islower() for i in stringa])-1 for stringa in self.sequences] # ok, this is likely to be unefficient but who cares...
                self.strip2align=_np.where([i!='.' and i!='*' and not i.islower() for i in self.sequences[0]])[0]
                #self.align2strip=_np.cumsum([i!='.' and i!='*' and not i.islower() for i in self.sequences[0]])-1
                self.align2strip=_np.ones(self.N_align,dtype=_np.int32)*-1
                self.align2strip[_np.where([i!='.' and i!='*' and not i.islower() for i in self.sequences[0]])[0]]=_np.arange(self.N)
                ###

                self.Z=_np.array([[self.letter2numer[aa] for aa in s] for s in self.stripped_seqs]) # this is a MxN _np.array with the stripped sequences as number # TODO: this is probably not so efficient but who cares

                self.q=_np.max(self.Z)+1 # for proteins is always 21, but we can keep it general in case somebody wants to use RNA


        def get_dict(self,stripped=False):
                """This creates a dictionary of names and sequences (analogous to the input of the constructor).
It may look useless. And it probably is."""
                if not stripped:
                        fasta_list=[ {'sequence':seq,'title':name} for name,seq in zip(self.names,self.sequences)]
                else:
                        fasta_list=[ {'sequence':seq,'title':name} for name,seq in zip(self.names,self.stripped_seqs)]
                return fasta_list

        def filter_gaps(self,limit):
                """These method filters sequences with more than a specific number of continuos gaps "-", defined by the user.
                It returns a new alignment object as an output!"""
                # TODO: all this is done in a very stupid way. but it was the first thing that came into my mind.
                assert limit>0, "'limit' should better be positive, don't you think?"
                lim_string='-'*limit
                fasta_list=self.get_dict() #creating a dictionary of the sequences
                print('Old number of sequences=%d'%self.M)
                for iseq in range(self.M-1,-1,-1):
                        if lim_string in self.stripped_seqs[iseq]:
                                fasta_list.pop(iseq)
                print('New number of sequences=%d'%len(fasta_list))
                new_alignment=Alignment(fasta_list)
                return new_alignment

        def cut_alignment(self,start_ndx,end_ndx):
                """Cuts an alignment object based on indexes referred to the stripped sequences.
                Returns a new alignment object as output!"""
                ### TODO: check that I'm not doing something too stupidly inefficient here...
                i1=self.strip2align[start_ndx]
                i2=self.strip2align[end_ndx]
                new_seqs=[''.join(list(seq)[i1:i2]) for seq in self.sequences]
                fasta_list=[ {'sequence':seq,'title':name} for name,seq in zip(self.names,new_seqs)]
                new_alignment=Alignment(fasta_list)
                return new_alignment

        def get_original_seq(self,iseq):
                """Returns a specific original sequence form the alignment"""
                seq=_np.array(list(self.sequences[iseq]))
                orig_seq=''.join(seq[self.orig2align[iseq]])
                return orig_seq
        
        def find_sequence(self,seq):
                """Looks through the original sequences to find a given sequence. Returns the corresponding index, name and aligned sequence."""
                results=[] ### can we have more than one sequences that match the query?
                for iseq in range(self.M):
                        orig_seq=self.get_original_seq(iseq)
                        # use uppercase to avoid ambiguity
                        if seq.upper()==orig_seq.upper():
                                results.append((iseq,self.names[iseq],self.sequences[iseq]))
                if len(results)==0:
                        return False
                if len(results)==1:
                        return results[0]
                if len(results)>1:
                        print("WARNING: there are %d sequences that match the one you provided. Maybe you have identical sequences in the alignment?"%len(results))
                        return results


        def find_name(self,name):
                """Looks through the alignment to find a sequence with a given name. Returns the corresponding index, name and aligned sequence."""
                results=[] ### can we have more than one sequences that match the query?
                for i,n in enumerate(self.names):
                        if name in n:
                                results.append((i,n,self.stripped_seqs[i]))
                if len(results)==0:
                        return False
                if len(results)==1:
                        return results[0]
                if len(results)>1:
                        print("WARNING: there are %d sequences that match the one you provided. Maybe you have identical sequences in the alignment?"%len(results))
                        return results

        def write_alignment(self,filename,mode='stripped'):
                fh=open(filename,'w')
                if mode=='stripped':
                        for n,s in zip(self.names,self.stripped_seqs):
                                fh.write(">%s\n%s\n"%(n,s))
                        fh.close()
                        return True
                if mode=='aligned':
                        for n,s in zip(self.names,self.sequences):
                                fh.write(">%s\n%s\n"%(n,s))
                        fh.close()
                        return True
                else:
                        print("ERROR: mode '%s' not supported"%mode)
                        fh.close()
                        return False
                

def read_alignment(inputfile,filter_limit=None,check_aminoacid=True,check_nucleicacid=False):
        """Reads an alignment from a .fasta format file and returns an Alignment object"""
        fasta_list=FASTA_parser(inputfile,check_aminoacid=check_aminoacid,check_nucleicacid=check_nucleicacid)
        alin=Alignment(fasta_list)
        if filter_limit!=None:
                alin=alin.filter_gaps(filter_limit) # TODO is this smart?
        return alin

### TODO: a function to write down an alignment in fasta format
