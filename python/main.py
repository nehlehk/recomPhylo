#======================================================================================================================
'''
Jan 2020
Bayesian phylogenetics inference code
This code is my first practice to understand the concept of phylogenetic reconstruction by using Bayesian.
Authors: Nehleh Kargarfard
'''
#======================================================================================================================
# import libraries
import sys
import math
from builtins import print
from symbol import import_as_name

import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import seaborn as sns
import numpy.linalg as la
import pprint
import dendropy





def get_base_frequencies(dna):
    base_freq =  {base: dna.count(base)/float(len(dna)) for base in 'ACGT'}

    return  base_freq
#===================================================================================================================
def avg_base_frequencies(alignment):
    n = len(alignment)
    a = c = g = t = 0
    for align in alignment:
       sum_bf =  get_base_frequencies(align)
       a += sum_bf['A']
       c += sum_bf['C']
       g += sum_bf['G']
       t += sum_bf['T']

    return [a/float(n) , c/float(n) , g/float(n) , t/float(n)]
#===================================================================================================================
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
def GTR(br_len,alignment ,pi = []*4 ,a = 1 , b = 1 , c = 1 , d = 1 , e = 1 , f = 1):
    n = alignment.sequence_size
    # n = len(alignment[0])
    mu = 0

    freq = np.zeros((4, 4))
    q = np.zeros((4, 4))
    sqrtPi =  np.zeros((4, 4))
    sqrtPiInv =  np.zeros((4, 4))
    exchang = np.zeros((4, 4))
    s = np.zeros((4, 4))
    fun = np.zeros(n)

    # pi = avg_base_frequencies(alignment)
    # pi = [0.25,0.25,0.25,0.25]
    freq = np.diag(pi)
    sqrtPi = np.diag(np.sqrt(pi))
    sqrtPiInv = np.diag(1.0/np.sqrt(pi))

    mu = 1 / (2 * (a * pi[0] * pi[1]) + (b * pi[0] * pi[2]) + (c * pi[0] * pi[3]) + (d * pi[1] * pi[2]) + (
                e * pi[1] * pi[3]) + (pi[2] * pi[3]))

    exchang[0][1] = exchang[1][0] = a
    exchang[0][2] = exchang[2][0] = b
    exchang[0][3] = exchang[3][0] = c
    exchang[1][2] = exchang[2][1] = d
    exchang[1][3] = exchang[3][1] = e
    exchang[2][3] = exchang[3][2] = f

    q = np.multiply(np.dot(exchang,freq), mu)


    for i in range(4):
        q[i][i] = -sum(q[i][0:4])


    s = np.dot(sqrtPi,np.dot(q,sqrtPiInv))

    eigval, eigvec = la.eig(s)
    eigvec_inv = la.inv(eigvec)

    left = np.dot(sqrtPi,eigvec)
    right = np.dot(eigvec_inv,sqrtPiInv)

    p = np.dot(left, np.dot(np.diag(np.exp(eigval * br_len)), right))

    k = 0
    dna1 = alignment[k]
    dna2 = alignment[k + 1]
    for index, base1 in enumerate(dna1):
       base2 = dna2[index]
       i = give_index(str(base1))
       j = give_index(str(base2))
       fun[index] =  pi[i] * p[i][j]


    ll = np.log(np.sum(fun))

    return ll
#===================================================================================================================
def max_likeGTR(alignment):

    initial_guess = 0.1,alignment,[0.28,0.22,0.20,0.3]
    result = spo.minimize(GTR , initial_guess ,method='nelder-mead' , options={'disp' : True}) #, options={'disp' : True}
    print("Result by using scipy:")
    print("theta = {} , MLE = {}".format(result.x, result.fun))
    return result
#=======================================================================================================================

alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/0_Research/PhD/Data/test.fasta"), schema="fasta" , )

max_likeGTR(alignment)

ll = GTR(alignment,0.1,[0.28,0.22,0.20,0.3],0.3,0.1,0.1,0.1,0.4,1)

print(ll)




#=======================================================================================================================
def ratio_matrix_computation(dna_list):
    i = 0
    seq_count = len(dna_list)
    seq_length = len(dna_list[0])
    ratio_matrix = {base :{base2 :0 for base2 in 'ACGT'} for base in 'ACGT'}

    # while (i < seq_count -1 ):
    dna1 = dna_list[i]
    dna2 = dna_list[i+1]
    for index, base1 in enumerate(dna1):
        base2 = dna2[index]
        ratio_matrix[base1][base2] += 1  #give count

    # ============================================================
    #dirty code :(
    for base1 in 'ACGT':
        for base2 in 'ACGT':
            ratio_matrix[base1][base2] = ratio_matrix[base1][base2] / seq_length # give ratio

        # i += 1

    return ratio_matrix
#===================================================================================================================
def  transition_probability_matrix_JC(t,beta):
    q = np.zeros((4, 4))
    p = np.zeros((4, 4))
    q[:][:] = beta
    q[0][0] = q[1][1] = q[2][2] = q[3][3] = -3*beta
    eigval, eigvec = la.eig(q)
    eigvec_inv = la.inv(eigvec)
    p = np.dot(eigvec,np.dot(np.diag(np.exp(eigval * t)), eigvec_inv))
    return p
#=======================================================================================================================

def plot_MLE_GTR():
    d = np.arange(0.0001, 1 , 0.01)
    ll_array = []
    for i in d:
        y_plot = max_likeGTR()
        ll_array.append(y_plot)

    plt.plot(d, ll_array, label='fitted model')
    plt.axvline(max_like(), color='r', ls='-.')
    plt.title(str(float(max_like()))+"  is the maximum likelihood estimate (MLE) of v")
    plt.ylabel('log(L)')
    plt.xlabel('v')
    #plt.legend(loc='lower right')
    plt.show()
#=======================================================================================================================


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
alignment = dendropy.DnaCharacterMatrix.get\
    (file=open("/home/nehleh/0_Research/PhD/Data/test.fasta"), schema="fasta" , )




# v = np.arange( 0 , 2, 0.01)
# y1 = []
# y2 = []
# d = []
# for i in v:
#     # p = transition_probability_matrix_JC(i,i)
#     transition_probability_matrix_GTR(alignment,i)
    # y1.append(p[0][0])
    # y2.append(p[0][1])
    # d.append(i * i)



# pprint.pprint(p)
# plt.plot(d, y1, label='P(AA)')
# plt.plot(d, y2, label='P(AC)')
# plt.ylabel('probablity')
# plt.xlabel('v')
# plt.axhline(p[0][0], color='r', ls='-.')
# plt.legend(loc='lower right')
# plt.show()






# max_like()