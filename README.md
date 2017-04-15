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

- check the math. I really don't know what is going on in the DI calculation. I'm trusting the matlab code.

- input format: the octave version reads "." characters (I think an alternative is "*"). We can easily use both, I guess.

- efficiency: analysis of the large sample dataset took 3 hours on shark.imp.fu-berlin.de with the octave code (but then, maybe matlab is faster). with pydca it takes a couple of minutes

- Would be nice to be able to use it for different purposes (Granata et al., etc) and maybe also have different options for the DCA implementation (mfDCA, plmDCA, etc.).

- dca class: output DI as array; method/lists to go from "absolute" indexes to "local" indexes for each seq of the input; something to sort the DI; some get_contact_map(n. of highest DI to consider); remember we may want to exclude neighbours on the chain (up to 4 residues) when defining the highest DIs, but reconsider them later.

6) some function to rapidly plot the contact map, or other stuff

7) we making this for jupyter?  Or we want also a "main" script, to execute it directly from the terminal? It shouldn't be hard to make it.

8) also, it may be useful to have some function/method to "filter" sequences with too many adiacent gaps (Faruck gave me an independent script to do it. But I would put everything together.) I think this should go in the class "Alignment"