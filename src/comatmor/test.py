# Just to test a few things
import numpy as np
from comatmor import stationRB
A = stationRB()
A.addMatrix('Stiffnessmatrix', np.array([0,1]), np.array([0,1]), np.array([1,1]), 'Diffusion', 1, np.array([1,3]), np.array([1,3]))
A.addRhs(np.array([1,1]))
A.constructRB()

#================= multi dim param
	                        # case of multi-dimensional parameter dependence, not working so far
                                else:
                                        dim = len(paramName)
                                        paramType = ParameterType({ paramName[i]: paramShape[i] for i in range(dim)})
                                        print paramType
                                        paramFunc = GenericParameterFunctional(lambda mu: np.diag([mu[paramName][0] for paramName in paramType]), parameter_type = paramType)

