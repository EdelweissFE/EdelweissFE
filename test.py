import numpy as np
from scipy.sparse import csr_matrix


from fe.linsolve.pardiso.pardiso import LinearSolver

solver = LinearSolver()

A = np.array([[1.0, 0.0, 0.],[0., 1., 0.], [0., 0., 1.]])
A = csr_matrix(A)

b =np.array( [ [1., 2, 3], [3., 2., 1]]).T

solver.analyseSystem(A)
x = solver.solve(A, b)

print(x)

