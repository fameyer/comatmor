#!../../../../virt/bin/python
# Script to start RB-computations with harddisc data access with matlab
# by Falk Meyer, 12.02.2015

"""startRB

Usage:
  startHeatRB.py

Arguments:
  None

Options:
  None
"""

from docopt import docopt

from comatmor.elliptic import ellipticRB

def startRB(args):
	"""
	DOC ME
	"""

	# create stationRB object
	rb = ellipticRB(inputmethod = 'disc')
			
	# compute rb solution
	rb.constructRB(10)

	# One could also save the RB object
	# rb.save()

	rb.compute()#training_set)

if __name__ == '__main__':
	# parse arguments
	args = docopt(__doc__)
	# run script
        startRB(args)
