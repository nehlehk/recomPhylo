# import libraries
import sys
import math
from builtins import print
from symbol import import_as_name

import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import numpy.linalg as la
import pprint
from scipy.optimize import Bounds
import dendropy


def give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3
#=======================================================================================================================
# The first parameter (i.e theta) contains model parameters i.e: in JC69  it just contains edge length v.
#=======================================================================================================================
def log_likelihood(theta, data = None):

    if data == None:
        x = 90
        n = 948
    else:
        x = data[0] # number of mismatch
        n = data[1] # length of alignment


    # ===================   for two sequences   =======================
    v = theta[0]
    # p0 =  0.25 + 0.75 * math.exp(-4 * v / 3)
    # p1 =  0.75 - 0.75 * math.exp(-4 * v / 3)
    p0 =  1/16 + 3/16 * math.exp(-4 * v / 3)
    p1 =  1/16 - 1/16 * math.exp(-4 * v / 3)
    y = (n-x)* np.log(p0) + x * np.log(p1)
    # negLL = -y
    # print("v = {} , y = {}  , negLL = {} ".format(v, y , negLL))  # for tracing
    return -y

#=======================================================================================================================
# Using scipy.optimize.minimize to find MLE
#=======================================================================================================================
def  JC(param):
    n = alignment.sequence_size
    q = np.zeros((4, 4))
    p = np.zeros((4, 4))
    fun = np.zeros(n)

    q[:][:] = 1/3
    q[0][0] = q[1][1] = q[2][2] = q[3][3] = -1
    # pprint.pprint(q)
    eigval, eigvec = la.eig(q)
    eigvec_inv = la.inv(eigvec)
    p = np.dot(eigvec,np.dot(np.diag(np.exp(eigval * param[0])), eigvec_inv))

    k = 0
    dna1 = alignment[k]
    dna2 = alignment[k + 1]
    for index, base1 in enumerate(dna1):
       base2 = dna2[index]
       i = give_index(str(base1))
       j = give_index(str(base2))
       fun[index] = pi[i] * p[i][j]

    pprint.pprint(p)

    ll = np.sum(np.log(fun))

    return -ll
#===================================================================================================================
def max_like():
    initial_guess = [0.5]
    # result = spo.minimize(log_likelihood , initial_guess ,method='nelder-mead' , options={'disp' : True}) #, options={'disp' : True}

    bounds = Bounds([0], [100])
    result = spo.minimize(JC, initial_guess, method='trust-constr',options={'disp': True} , bounds =bounds)
    print("Result by using scipy:")
    print("d = {} , MLE = {}".format(result.x, -result.fun))
    return result.x

#=======================================================================================================================
# Using simple loop to find MLE
#=======================================================================================================================
def max_like_manual():
    v = np.arange( 0.01, 1, 0.001)
    ll_array = []
    for i in v:
        y_plot = log_likelihood([i])
        ll_array.append(y_plot)

    MLE = np.min(ll_array)
    result = np.where(ll_array == np.amin(ll_array))
    print("Result by using simple loop:")
    print("theta:", float(v[result]), "MLE:", MLE)
    return float(v[result])
#=======================================================================================================================
# Plot Maximum Likelihood Estimation
#=======================================================================================================================
def plot_MLE():
    d = np.arange(0.0001, max_like() + 0.2, 0.005)
    ll_array = []
    for i in d:
        y_plot = -log_likelihood([i])
        ll_array.append(y_plot)

    plt.plot(d, ll_array, label='fitted model')
    plt.axvline(max_like(), color='r', ls='-.')
    plt.title(str(float(max_like()))+"  is the maximum likelihood estimate (MLE) of v")
    plt.ylabel('log(L)')
    plt.xlabel('v')
    #plt.legend(loc='lower right')
    plt.show()

#=======================================================================================================================

# max_like_manual()

pi = [0.25]*4

alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/0_Research/PhD/Data/test.fasta"), schema="fasta")

max_like()



# plot_MLE()