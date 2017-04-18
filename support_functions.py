import re
import datetime
import numpy as np

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
                self.numer2letter=list(self.letter2numer.keys())[:-1] # NB: 0 will always be backmapped to "-"
                
                self.M=len(fasta_list) # this is the number of sequences
                
                self.sequences=[seq['sequence'] for seq in fasta_list] # sequences
                self.names=[seq['title'] for seq in fasta_list] # names of sequences can be useful in the prostprocessing

                self.orig2align=[np.where([i!='-' or i!='.' for i in stringa])[0] for stringa in self.sequences] # this takes x3 times than stripping the seq.
                self.align2orig=[np.cumsum([i!='-' or i!='.' for i in stringa])-1 for stringa in self.sequences] # ok, this is likely to be unefficient but who cares...
                
                self.__strip()

                
        def __strip(self):
                """Fuck this. I'm going to be a stripper"""
                this_name='stripper'
                self.stripped_seqs=[re.sub("[a-z]|\.|\*",'',s) for s in self.sequences] # this should remove lowercase letters and "."  and "*" ### TODO: do we really need to do search for . and * for all sequences??
                ### TODO: these two lists of lists are probably redundant. Since the gaps are always in the same positions only one list for all the sequences is enough...
                self.strip2align=[np.where([i!='.' or i!='*' or i.islower() for i in stringa])[0] for stringa in self.sequences] # this takes x3 times than stripping the seq.
                self.align2strip=[np.cumsum([i!='.' or i!='*' or i.islower() for i in stringa])-1 for stringa in self.sequences] # ok, this is likely to be unefficient but who cares...

                if len(set([len(s) for s in self.stripped_seqs]))>1:
                        raise_error(this_name,"ERROR: stripped sequences have different lengths!")
                self.N=len(self.stripped_seqs[0]) # this is the lenght of each stripped sequences

                self.Z=np.array([[self.letter2numer[aa] for aa in s] for s in self.stripped_seqs]) # this is a MxN np.array with the stripped sequences as numbers

                self.q=np.max(self.Z)+1 # for proteins is always 21, but we can keep it general in case somebody wants to use RNA
