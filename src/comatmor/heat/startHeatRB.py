#!../../../../virt/bin/python
# Script to start RB-computations with harddisc data access with matlab
# by Falk Meyer, 12.02.2015

"""startHeatRB

Usage:
  startHeatRB.py [--endtime=TIME] [--steps=STEPS] [--samples=SAMPLES]

Arguments:
  None

Options:
  --endtime=TIME	Use TIME as end-time for the instationary problem

  --steps=STEPS		STEPS to be performed during timestepping.

  --samples=SAMPLES	Number of sample parameters, the reduced basis shall be constructed for

"""

from docopt import docopt

from comatmor.heat import instationHeatRB

def startHeatRB(args):
	"""
	Script to start RB computation for given time and time steps, intended to be called by matlab.
	"""
	T = float(args['--endtime'] or 1)
	step_number = int(args['--steps'] or 10)
	num_samples = int(args['--samples'] or 10)
	# create stationRB object
	rb = instationHeatRB(inputmethod = 'disc')
			
	# construct RB basis for num_samples in parameterspace, n steps and end-time T
	rb.constructRB(num_samples,T,step_number)

	# Compute solutions for given training_set	
	rb.compute()

if __name__ == '__main__':
	# parse arguments
	args = docopt(__doc__)
	# run script
        startHeatRB(args)
