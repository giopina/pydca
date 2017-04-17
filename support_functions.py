import re
import datetime
import numpy as np

def write_log(name, text):
        log_filename = name + '.log'
        log_file = open(log_filename, 'w')
        log_file.write(text)
        log_file.close()
        

def header(name):
	return "[" + name + " " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "] "


def raise_error(name, text):
	indent = " "*len(header(name))
	lines = text.split('\n')
	logmsg = header(name) + lines[0] + '\n'
	for line in lines[1:]:
		logmsg += indent + line + '\n'
	write_log(name, logmsg)
	raise NameError(logmsg)


def print_log(name, text, quiet=False):
	indent = " "*len(header(name))
	lines = text.split('\n')
	logmsg = header(name) + lines[0]
	for line in lines[1:]:
		logmsg += '\n' + indent + line
	if not quiet:
		print(logmsg)
	return logmsg + '\n'


def FASTA_parser(fasta_filename, to_upper=False, check_aminoacid=False, check_nucleicacid=False):
	this_name = "FASTA_parser"

	if check_aminoacid and check_nucleicacid:
		raise_error(this_name, "ERROR: only one out of 'check_aminoacid' and 'check_nucleicacid' options can be active")

	fasta_list = []
	fasta_file = open(fasta_filename, 'r')
	text = fasta_file.read().split('\n')
	fasta_file.close()

	there_is_title_line = False
	nl = 0
	no_titles = 0
	for t_line in text:
		nl += 1
		line = t_line.strip()
		if not line:
			there_is_title_line = False
			continue
		if not (line[0] == '>' or re.match('^[A-Za-z\-\*\\\\]+$', line)):
			raise_error(this_name, "ERROR: Line {0} does not respect the FASTA format".format(nl))
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
		elif re.match('^[A-Za-z\-\*\\\\]+$', line):
			if not fasta_list or (fasta_list[-1]['sequence'] == '' and not there_is_title_line):
				raise_error(this_name, "ERROR: there is a string of text not preceded by a title at line {0}".format(nl))
			if to_upper:
				fasta_list[-1]['sequence'] += line.strip().replace('\\', '').upper()
			else:
				fasta_list[-1]['sequence'] += line.strip().replace('\\', '')
			there_is_title_line = False

	allowed_characters = None
	if check_aminoacid:
		allowed_characters = 'ABCDEFGHIKLMNPQRSTUVWYZXabcdefghiklmnpqrstuvwyzx*-'
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
                TODO: add a method to go back to the original indexing of each sequence from the stripped indexes
                """
                self.letter2numer={\
                # full AA alphabet
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
                }
                self.numer2letter=self.letter2numer.keys()[:-1] # NB: 0 will always be backmapped to "-"
                
                self.M=len(fasta_list) # this is the number of sequences
                
                self.sequences=[seq['sequence'] for seq in fasta_list] # sequences
                self.names=[seq['title'] for seq in fasta_list] # names of sequences can be useful in the prostprocessing
                self.__strip()

                
        def __strip(self):
                """Fuck this. I'm going to be a stripper"""
                this_name='stripper'                
                self.stripped_seqs=[re.sub("[a-z]|\.|\*",'',s) for s in self.sequences] # this should remove lowercase letters and "."  and "*"
                #stripped_seqs=[np.cumsu[i!='.' and i!='-' for i in stringa] for stringa in self.sequences] ### this is slower by a factor x2
                self.strip2align=[np.where([i!='.' and i!='-' for i in stringa]) for stringa in self.sequences] # this takes x3 times than stripping the seq.
                self.align2strip=[np.cumsum([i!='.' and i!='-' for i in stringa])-1 for stringa in self.sequences] # ok, this is likely to be unefficient but who cares...

                if len(set([len(s) for s in self.stripped_seqs]))>1:
                        raise_error(this_name,"ERROR: stripped sequences have different lengths!")
                self.N=len(self.stripped_seqs[0]) # this is the lenght of each stripped sequences

                self.Z=np.array([[self.letter2numer[aa] for aa in s] for s in self.stripped_seqs]) # this is a MxN np.array with the stripped sequences as numbers

                self.q=np.max(self.Z)+1 # for proteins is always 21, but we can keep it general in case somebody wants to use RNA
