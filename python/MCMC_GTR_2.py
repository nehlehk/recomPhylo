# import libraries
import sys
import math
from builtins import print
from symbol import import_as_name

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import dendropy
import numpy.linalg as la


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
def GTR(params):
    n = alignment.sequence_size
    mu = 0

    freq = np.zeros((4, 4))
    q = np.zeros((4, 4))
    sqrtPi =  np.zeros((4, 4))
    sqrtPiInv =  np.zeros((4, 4))
    exchang = np.zeros((4, 4))
    s = np.zeros((4, 4))
    fun = np.zeros(n)
    a,b,c,d,e  = params[1:6]
    f = 1

    pi = params[6:10]

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

    p = np.dot(left, np.dot(np.diag(np.exp(eigval * params[0])), right))

    k = 0
    dna1 = alignment[k]
    dna2 = alignment[k + 1]
    for index, base1 in enumerate(dna1):
       base2 = dna2[index]
       i = give_index(str(base1))
       j = give_index(str(base2))
       fun[index] = pi[i] * p[i][j]

    ll = np.sum(np.log(fun))

    return ll
#======================================================================================================================
#The tranistion model defines how to move from current to new
def transition_model(theta,type,width):
    proposal = np.zeros(len(theta))
    if type == 'sw_norm':
        for index,param in enumerate(theta):
            proposal[index] = np.random.normal(param, width)
            if (proposal[index] < 0):
                proposal[index] =  -param - proposal[index]  # reflection in normal dist --- the excess is reflected back into the interval
    # print(proposal)
    return proposal
#=======================================================================================================================
#Define prior
def prior(theta):
    return 1
#=======================================================================================================================
#Defines whether to accept or reject the new sample
def acceptance(current_pos, new_pos):
    if new_pos > current_pos:
        return True
    else:
        accept = np.random.uniform(0,1)
        # Since we did a log likelihood, we need to exponentiate in order to compare to the random number less likely
        # new_pos are less likely to be accepted
        return  (accept < (np.exp(new_pos - current_pos)))
    # for index, param in enumerate(new_pos):
    #     if new_pos[index] > current_pos[index]:
    #         res[index] = 1
    #     else:
    #         accept = np.random.uniform(0,1)
    #         if (accept < (np.exp(new_pos[index] - current_pos[index]))):
    #             res[index] = 1
    #
    # return sum(res) == len(current_pos)
# =======================================================================================================================
def metropolis_hastings(likelihood_computer,prior, transition_model, param_init,iterations,data,acceptance_rule):
    # likelihood_computer(theta,data): returns the likelihood that these parameters generated the data
    # transition_model(theta): a function that draws a sample from a symmetric distribution and returns it
    # param_init: a starting sample
    # iterations: number of accepted to generated
    # data: the data that we wish to model
    # acceptance_rule(current_pos, new_pos): decides whether to accept or reject the new sample

    temp = len(param_init)-2
    theta = param_init[:temp]
    proposal_type = param_init[temp]
    proposal_width = param_init[temp+1]
    accepted = []
    rejected = []
    final_log_acc = []
    final_log_rej = []
    all_theta = []
    for i in range(iterations):
        new_theta = transition_model(theta,proposal_type,proposal_width)
        lik_theta = likelihood_computer(theta)
        lik_new_theta = likelihood_computer(new_theta)
        if (acceptance(lik_theta+np.log(prior(theta)), lik_new_theta+np.log(prior(new_theta)))):
            theta = new_theta
            accepted.append(new_theta)
            if not(i % 1):
                final_log_acc.append(new_theta)
        else:
            rejected.append(new_theta)
            if not(i % 1):
                final_log_rej.append(new_theta)
        all_theta.append(theta)


    return np.array(final_log_acc) , np.array(final_log_rej) , np.array(accepted) , np.array(rejected) , np.array(all_theta)
# ======================================================================================================================
# pi = [0.25]*4
iterations = 5000

alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/0_Research/PhD/Data/test.fasta"), schema="fasta")

final_log_acc,final_log_rej,accepted,rejected,all_theta  =  \
    metropolis_hastings(GTR,prior,transition_model,[0.2,1,1,1,1,1,1/4,1/4,1/4,1/4,'sw_norm',0.01],iterations,alignment,acceptance)

print("Final_Acc",final_log_acc.shape)
print("Final_rej",final_log_rej.shape)



avg = np.zeros(final_log_acc.shape[1])
for item in range(final_log_acc.shape[1]):
    avg[item] = np.mean(final_log_acc[:,item])



print(avg)



# print("edge_len = {}, a={} , b={} , c ={} , d ={} , e={}".format(avg[0] ,avg[1],avg[2],avg[3],avg[4],avg[5] ))

# print("P_jump = ", round(accepted.shape[0]/iterations * 100 , 2) , '%' )







# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(2,1,1)
# ax.plot( rejected[:50], 'rx', label='Rejected',alpha=0.9)
# ax.plot( accepted[:50], 'b.', label='Accepted',alpha=0.9)
# ax.set_xlabel("Iteration")
# ax.set_ylabel("v")
# ax.set_title("Figure 2: MCMC sampling for v with Metropolis-Hastings. First 500 samples are shown.")
# ax.grid()
# ax.legend()
#
# ax2 = fig.add_subplot(2,1,2)
# to_show=final_log_acc.shape[0]
# ax2.plot( final_log_rej[:to_show], 'rx', label='Rejected',alpha=0.5)
# ax2.plot( final_log_acc[:to_show], 'b.', label='Accepted',alpha=0.5)
# ax2.set_xlabel("Iteration")
# ax2.set_ylabel("v")
# ax2.set_title("Figure 3: MCMC sampling for v with Metropolis-Hastings. All samples are shown.")
# ax2.grid()
# ax2.legend(loc="best")
#
# fig.tight_layout()
#
#
# temp=int(0.25*all_theta.shape[0])
# print("Mean_whole_theta after burn-in: ",np.mean(all_theta[temp:]))
# print("SD_whole_theta after burn-in: ",np.std(all_theta[temp:]))
#
# # consider the initial 25% of the values of v to be "burn-in", so we drop them.
# show=int(0.25*final_log_acc.shape[0])
#
# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(1,1,1)
# ax.plot(final_log_acc[show:])
# ax.set_title("Figure 4: Trace for v")
# ax.set_ylabel("v")
# ax.set_xlabel("Iteration")
#
#
#
#
#
# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(1,1,1)
# # ax.hist(accepted,bins=300,label="Predicted density")
# ax.set_xlabel("$\Theta$")
# ax.set_ylabel("Frequency")
# ax.set_title("Figure 5: posterior densities for sequence distance $\Theta$ under the JC69 model")
# # ax.legend()
# ax.grid("off")
# ax = sns.distplot(accepted,bins=300,hist= True)
#
# plt.show()