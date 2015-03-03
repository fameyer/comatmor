#!../../../../virt/bin/python
# Script to start RB-computations with harddisc data access with matlab
# by Falk Meyer, 12.02.2015

"""startRB

Usage:
  startHeatRB.py [--samples=SAMPLES]

Arguments:
  None

Options:
  --samples=SAMPLES     Number of sample parameters, the reduced basis shall be constructed for

"""

from docopt import docopt

from comatmor.elliptic import ellipticRB

def startRB(args):
	"""
	DOC ME
	"""
	num_samples = int(args['--samples'] or 10)

	# create stationRB object
	rb = ellipticRB(inputmethod = 'disc')
			
	# compute rb solution
	rb.constructRB(num_samples)

	rb.compute()

if __name__ == '__main__':
	# parse arguments
	args = docopt(__doc__)
	# run script
        startRB(args)
