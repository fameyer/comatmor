#!../../../virt/bin/python
# Script to start RB-computations with harddisc data access with matlab
# by Falk Meyer, 12.02.2015

#import os
import comatmor

def startRB():
	"""
	DOC ME
	"""

	# create stationRB object
	rb = comatmor.stationRB(inputmethod = 'disc')
	
	# add system
	rb.addMatrix()
	
	# add rhs
	rb.addRhs()
	
	# compute rb solution
	rb.constructRB()

	# One could also save the RB object
	# rb.save()

	# compute several solutions in rb setting and write to harddisc
	training_set =[1.5, 2.0]
	
	rb.compute(training_set = training_set)

if __name__ == '__main__':
        startRB()

