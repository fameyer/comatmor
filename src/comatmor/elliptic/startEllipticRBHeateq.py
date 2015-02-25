#!../../../../virt/bin/python
# Script to start RB-computations with harddisc data access with matlab
# by Falk Meyer, 12.02.2015

from comatmor.elliptic import ellipticRB

def startRB():
	"""
	DOC ME
	"""

	# create stationRB object
	rb = ellipticRB(inputmethod = 'disc')
			
	# add system
	rb.addMatrix()
	
	# compute rb solution
	rb.constructRB(10)

	# One could also save the RB object
	# rb.save()

	# compute several solutions in rb setting and write to harddisc
	training_set =[(i,j) for i in range(1,50,5) for j in range(1,50,5)]
	
	rb.compute(training_set = [(1.0,40.0)])#training_set)

if __name__ == '__main__':
        startRB()

