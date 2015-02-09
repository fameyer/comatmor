# Main file to test COMSOL-MATLAB-PYMOR interface
# by Falk Meyer, 09.02.2015

# Standard import
import numpy as np

# pymor includes
from pymor.operators.numpy import NumpyMatrixOperator

from pymor.discretizations.basic import StationaryDiscretization

# Input routines
from input import getMatrix



def main():
	"""
	Doc me
	"""
	# get input matrices
	MatDic = getMatrix()
	print MatDic

	# for any given matrix create a NumpyMatrixOperator in a dict
	OperatorDic = {} 
	for key in MatDic:
		print 'Transform '+key+' into an operator...' 
		OperatorDic[key] = NumpyMatrixOperator(MatDic[key])
	print OperatorDic

	# CREATE LINCOMBOPERATOR

	# make it for a stationary case for example
	rhs = NumpyMatrixOperator(np.ones(1217))	
	for key in OperatorDic:
		stat = StationaryDiscretization(operator=OperatorDic[key],rhs=rhs)
		print stat.solve()
	 
	
if __name__ == '__main__':
	main()
