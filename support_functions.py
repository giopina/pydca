import re
import datetime

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
