#!../../../../virt/bin/python
# CHANGE ABOVE PYTHON REFERENCE TO YOUR lOCAL ONE

#    Copyright (C) 2015 Falk Meyer
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#===========================================================================

# Script to start RB-computations with harddisc data access with matlab
# by Falk Meyer, 12.02.2015

"""startHeatRB

Usage:
  startHeatRB.py [--endtime=TIME] [--steps=STEPS] [--samples=SAMPLES] [--max_extensions=MAX] [--target_error=ERROR]

Arguments:
  None

Options:
  --endtime=TIME	Use TIME as end-time for the instationary problem

  --steps=STEPS		STEPS to be performed during timestepping.

  --samples=SAMPLES	Number of sample parameters, the reduced basis shall be constructed for

  --max_extensions=MAX   MAX number of basis extensions in the reduced basis construction

  --target_error=ERROR  Target error to reach during reduced basis construction

"""

from docopt import docopt

from comatmor.IRT import IRTRB

def startIRTRB(args):
	"""
	Script to start RB computation for given time and time steps, intended to be called by matlab.
	"""
	T = float(args['--endtime'] or 5)
	step_number = int(args['--steps'] or 20)
	num_samples = int(args['--samples'] or 4)
	max_extensions = int(args['--max_extensions'] or 10)
	target_error = float(args['--target_error'] or 1e-10)
	# create stationRB object
	rb = IRTRB(inputmethod = 'disc')
			
	# construct RB basis for num_samples in parameterspace, n steps and end-time T
	rb.constructRB(num_samples,T,step_number, max_extensions, target_error)

	# Compute solutions for given training_set	
	rb.compute()

if __name__ == '__main__':
	# parse arguments
	args = docopt(__doc__)
	# run script
        startIRTRB(args)
