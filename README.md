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

- ALIGNMENT: the octave version reads "." characters (I think an alternative is "*"), we can easily use both, I guess;
method/lists to go from "absolute" indexes to "local" indexes for each seq of the input (this is done but should be tested maybe);
method to "filter" sequences with too many adiacent gaps (Faruck gave me an independent script to do it);

- DCA: output DI as array (actually is probably more important the list of the highest-DI pairs);
some get_contact_map(n. of highest DI to consider);
remember we may want to exclude neighbours on the chain (up to 4 residues) when defining the highest DIs, but may want to reconsider them later (TODO).

- efficiency: analysis of the large sample dataset took 3 hours on shark.imp.fu-berlin.de with the octave code (but then, maybe matlab is faster). with pydca it takes a couple of minutes

- memory: storing all the matrices (Pij, Z, takes a lot of memory (2GB for 478 residues of stripped_286). Does it really make sense to store them after the DI is computed??

- Would be nice to be able to use it for different purposes (Granata et al., etc) and maybe also have different options for the DCA implementation (mfDCA, plmDCA, etc.)

- PLOT: some function to rapidly plot the contact map, or other stuff


