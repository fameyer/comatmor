#!../../../../virt/bin/python
# Script to start RB-computations with harddisc data access with matlab
# by Falk Meyer, 12.02.2015

"""startHeatRB

Usage:
  startHeatRB.py [--endtime=TIME] [--steps=STEPS]

Arguments:
  None

Options:
  --endtime=TIME	Use TIME as end-time for the instationary problem

  --steps=STEPS		STEPS to be performed during timestepping.
"""

from docopt import docopt

from comatmor.heat import instationHeatRB

def startHeatRB(args):
	"""
	Script to start RB computation for given time and time steps, intended to be called by matlab.
	"""
	T = float(args['--endtime'] or 1)
	step_number = int(args['--steps'] or 10)
	# create stationRB object
	rb = instationHeatRB(inputmethod = 'disc')
			
	# compute rb solution to with n steps and end-time T
	rb.constructRB(T,step_number)

	#print rb.getRB()

	# One could also save the RB object
	# rb.save()

	# compute several solutions in rb setting and write to harddisc
	#training_set =[(i,j) for i in range(1,50,5) for j in range(1,50,5)]
	
	rb.compute(training_set = [(1.0,40.0)])

if __name__ == '__main__':
	# parse arguments
	args = docopt(__doc__)
	# run script
        startHeatRB(args)
