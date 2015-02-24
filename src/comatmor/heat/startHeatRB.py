#!../../../../virt/bin/python
# Script to start RB-computations with harddisc data access with matlab
# by Falk Meyer, 12.02.2015

from comatmor.heat import instationHeatRB

def startRB():
	"""
	DOC ME
	"""

	# create stationRB object
	rb = instationHeatRB(inputmethod = 'disc')
			
	# add system
	rb.addMatrix()
	
	# add mass matrix
	rb.getMass()

	# add rhs, deprecated
	#rb.addRhs()

	# add initial solution
	rb.getU0()
	
	# compute rb solution to with n steps and end-time T
	T = 1 
	n = 10
	rb.constructRB(T,n)

	#print rb.getRB()

	# One could also save the RB object
	# rb.save()

	# compute several solutions in rb setting and write to harddisc
	training_set =[(i,j) for i in range(1,50,5) for j in range(1,50,5)]
	
	rb.compute(training_set = [(1.0,1.0)])#training_set)

if __name__ == '__main__':
        startRB()

