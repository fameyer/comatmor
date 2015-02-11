# Just to test a few things
import numpy as np
from comatmor import stationRB
A = stationRB()
A.addMatrix('Stiffnessmatrix', np.array([0,1]), np.array([0,1]), np.array([1,1]), 'Diffusion', 1, np.array([1,3]), np.array([1,3]))
A.addRhs(np.array([1,1]))
A.constructRB()
