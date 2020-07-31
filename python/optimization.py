import numpy as np
from scipy import optimize
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint


def rosen(x):
    return 100*(x[1]-x[0]**2)**2 + (1-x[0])**2


def cons_f(x):
    return [x[0]**2 + x[1], x[0]**2 - x[1]]


def cons_J(x):
    return [[2*x[0], 1], [2*x[0], -1]]


def cons_H(x, v):
    return v[0]*np.array([[2, 0], [0, 0]]) + v[1]*np.array([[2, 0], [0, 0]])


nonlinear_constraint = NonlinearConstraint(cons_f, -np.inf, 1, jac=cons_J, hess=cons_H)


bounds = Bounds([0, -0.5], [1.0, 2.0])
linear_constraint = LinearConstraint([[1, 2], [2, 1]], [-np.inf, 1], [1, 1])

x0 = np.array([0.5, 0])

res = optimize.minimize(rosen, x0, method='trust-constr',
                        constraints=[linear_constraint, nonlinear_constraint],
                        options={'verbose': 1}, bounds=bounds)

print(res.x)
