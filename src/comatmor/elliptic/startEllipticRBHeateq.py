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

"""startRB

Usage:
  startHeatRB.py [--samples=SAMPLES] [--max_extensions=MAXE] [--target_error=ERROR]

Arguments:
  None

Options:
  --samples=SAMPLES    	Number of sample parameters, the reduced basis shall be constructed for

  --max_extensions=MAX	MAX number of basis extensions in the reduced basis construction

  --target_error=ERROR  Target error to reach during reduced basis construction

"""

from docopt import docopt

from comatmor.elliptic import ellipticRB

def startRB(args):
	"""
	DOC ME
	"""
	num_samples = int(args['--samples'] or 10)
	max_extensions = int(args['--max_extensions'] or 30)
	target_error = float(args['--target_error'] or 1e-10)

	# create stationRB object
	rb = ellipticRB()
			
	# compute rb solution
	rb.constructRB(num_samples, max_extensions, target_error)

	rb.compute()

if __name__ == '__main__':
	# parse arguments
	args = docopt(__doc__)
	# run script
        startRB(args)
