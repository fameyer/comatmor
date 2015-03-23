#!../../../../virt/bin/python

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
