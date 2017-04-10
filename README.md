# pydca
Python implementation of mean-field direct coupling analysis
see
http://dca.rice.edu/portal/dca/
for references and original code in Matlab

## subfolders
##### Octave
Octave working version (tested)
##### test
Scripts to test against one dataset (reference results produced with the original matlab code)

## TODO:
1) read fasta: create class "alignment", that mimics what is used in dca.m
(fastaread() returns a list of "alignment" objects that have an attribute .Sequence. Maybe it's not necessary to build a class for it. Probably is good for the future)

2) input format: the octave version reads "." characters (I think an alternative is "*"). We can easily use both, I guess.

3) efficiency: it takes something like hours to run on FU workstations. Maybe is the algorithm itself that is computationally expensive, but let's be careful not to make any obvious mistake that consumes time. Let's also try to understand what's the time consuming part.

4) portability: let's create a something like a class "dca" (what about an API? actually, I'm not sure how to do it correcly and clean so I'd go for something rough but functional).
Would be nice to be able to use it for different purposes (Granata et al., etc) and maybe also have different options for the DCA implementation (mfDCA, plmDCA, etc.).

5) dca class: must contain: input alignment; output DI; method/lists to go from "absolute" indexes to "local" indexes for each seq of the input; something to sort the DI; some get_contact_map(n. of highest DI to consider); remember we may want to exclude neighbours on the chain (up to 4 residues) when defining the highest DIs, but reconsider them later.

6) some function to rapidly plot the contact map

7) we making this for jupyter?  Or we want also a "main" script, to execute it directly from the terminal? It shouldn't be hard to make it.

8) also, it may be useful to have some function/method to "filter" sequences with too many adiacent gaps (Faruck gave me an independent script to do it. But I would put everything together.)