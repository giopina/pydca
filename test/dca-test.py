import pydca.dca as dca
import sys
infile=sys.argv[1]
outfile=sys.argv[2]
dca_obj=dca.compute_dca(infile,compute_MI=True)
dca_obj.print_results(outfile)
