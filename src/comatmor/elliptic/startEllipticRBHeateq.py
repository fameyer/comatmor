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
