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

- ALIGNMENT: method/lists to go from "absolute" indexes to "local" indexes for each seq of the input (this is done but should be tested maybe. GP: now pretty much tested);

- DCA: output DI as array (actually is probably more important the list of the highest-DI pairs);
       some get_contact_map(n. of highest DI to consider);
       remember we may want to exclude neighbours on the chain (up to 4 residues) when defining the highest DIs, but may want to reconsider them later (TODO).

- memory: storing all the matrices (Pij, C, take a lot of memory (2GB for 478 residues of stripped_286). maybe exploiting symmetries one can save 50% memory (is it worth it?)

- Would be nice to be able to use it for different purposes (Granata et al., etc). What can we do? Implement the eSPECTRUS method, or is it too much?

- TODO: different options for the DCA implementation (mfDCA, plmDCA, etc.)

- PLOT: some function to rapidly plot the contact map, or other stuff


