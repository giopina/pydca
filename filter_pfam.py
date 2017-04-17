######################################################################################
## This program filters and cleans Pfam domain sequences taking of the dots "." and ##
## lower case residues, in order to reduce file size. These positions are ignored   ##
## during DCA analysis. Also this program filters bad sequences with more than a    ##
## specific number of continuous gaps "-", defined by the user.   		    ##	
######################################################################################

import linecache
import textwrap

data = raw_input("Pfam file name: ")
size = open(data,"r")
limit = int(raw_input("Maximum number of continous gaps per sequence: "))

###################################################################################
lim =""
for k in range(0,limit):
	lim+="-"
#print lim
###################################################################################


i=1

l = len(size.readlines())


output = open(data+"_filtered"+str(limit),"w")

####################################################################################

nseq = 0
excluded =0

####################################################################################

while i < l:
	sequence = ""
	n = linecache.getline(data, i)	
	counter = 0	
	if n[0] == ">":
		name = n
		next =  linecache.getline(data, i+1)
		try:
			while next[0] != ">":
				sequence=sequence+next
				i+=1
				next =  linecache.getline(data, i+1)
			nseq+=1
		except IndexError:
			pass
	i+=1
	x=""
	for j in range(0,len(sequence)):
		if sequence[j]!="." and sequence[j]!="\n" and sequence[j].islower()==False:
			x+=sequence[j]
	if len(x.split(lim)) == 1:
		output.write(name+x+"\n")
	else:
		excluded+=1

output.close()


print "\nOriginal number of sequences: "+str(nseq)+"\n"
percentage = float(excluded)/nseq*100

print "\nNumber(%) of sequences excluded: "+str(excluded)+" ("+str("{0:.2f}".format(percentage))+"%)\n"

print "\nFile saved as: "+data+"_filtered"+str(limit)+"\n"



