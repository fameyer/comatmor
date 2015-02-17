#!../../../virt/bin/python
# Script to start RB-computations with harddisc data access with matlab
# by Falk Meyer, 12.02.2015

import os
import comatmor

def startRB():
	"""
	DOC ME
	"""
	# Call comsol server matlab with extraction of everything necessary
	print 'Starting comsol-matlab server...'
 	os.system('comsol server matlab < basic.m')
	# commands to connect local comsol to server (when model is not present)	
	# commands to run basic.m
	# perhaps all in one

	# create stationRB object
	rb = comatmor.stationRB(inputmethod = 'disc')
	
	# add system
	rb.addMatrix()
	
	# add rhs
	rb.addRhs()
	
	# compute rb solution
	rb.constructRB()

	# compute several solutions in rb setting and write to harddisc
	training_set =[1.5, 2.0]
	
	rb.compute(training_set = training_set)

if __name__ == '__main__':
        startRB()

