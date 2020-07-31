# import libraries
import sys
import math
from builtins import print
from symbol import import_as_name

import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import seaborn as sns
import pprint



def log_likelihood(theta, data = None, model = 'JC69'):

    if data == None:
        x = 90
        n = 948
    else:
        x = data[0] # number of mismatch
        n = data[1] # length of alignment

    if model == 'JC69':
        # ===================   for two sequences   =======================
        v = theta
        # p0 =  0.25 + 0.75 * math.exp(-4 * v / 3)
        # p1 =  0.75 - 0.75 * math.exp(-4 * v / 3)
        p0 = 1 / 16 + 3 / 16 * math.exp(-4 * v / 3)
        p1 = 1 / 16 - 1 / 16 * math.exp(-4 * v / 3)
        y = (n-x)* np.log(p0) + x * np.log(p1)
        # negLL = -y
        # print("v = {} , y = {}  , negLL = {} ".format(v, y , negLL))  # for tracing
        return -y
#======================================================================================================================
#The tranistion model defines how to move from current to new
def transition_model(theta,type,width):
    if type == 'sw_norm':
        proposal = np.random.normal(theta, width)
        if (proposal < 0):
            proposal =  -theta - proposal  # reflection in normal dist --- the excess is reflected back into the interval
    elif type == 'sw_uni':
        proposal = np.random.uniform(theta - width/2 , theta + width/2)
        if (proposal < 0):
            proposal = - proposal  # reflection in uniform dist
    elif type == 'pscale':
        c = np.exp(width * (np.random.uniform(0,1) -1)/2 )
        proposal = theta * c


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
# =======================================================================================================================
def metropolis_hastings(likelihood_computer,prior, transition_model, param_init,iterations,data,acceptance_rule):
    # likelihood_computer(theta,data): returns the likelihood that these parameters generated the data
    # transition_model(theta): a function that draws a sample from a symmetric distribution and returns it
    # param_init: a starting sample
    # iterations: number of accepted to generated
    # data: the data that we wish to model
    # acceptance_rule(current_pos, new_pos): decides whether to accept or reject the new sample


    theta = param_init[0]
    proposal_type = param_init[1]
    proposal_width = param_init[2]
    accepted = []
    rejected = []
    final_log_acc = []
    final_log_rej = []
    all_theta = []
    for i in range(iterations):
        new_theta = transition_model(theta,proposal_type,proposal_width)
        lik_theta = likelihood_computer(theta,data)
        lik_new_theta = likelihood_computer(new_theta,data)
        if (acceptance(lik_theta+np.log(prior(theta)), lik_new_theta+np.log(prior(new_theta)))):
            theta = new_theta
            accepted.append(new_theta)
            if not(i % 100):
                final_log_acc.append(new_theta)
        else:
            rejected.append(new_theta)
            if not(i % 100):
                final_log_rej.append(new_theta)
        all_theta.append(theta)


    return np.array(final_log_acc) , np.array(final_log_rej) , np.array(accepted) , np.array(rejected) , np.array(all_theta)
# ======================================================================================================================

iterations = 500000

# transition_model(0.5,'pscale',0.1)

final_log_acc,final_log_rej,accepted,rejected,all_theta  =  \
    metropolis_hastings(log_likelihood,prior,transition_model,[0.2,'sw_norm',0.1],iterations,[90,948],acceptance)



print("P_jump = ", round(accepted.shape[0]/iterations * 100 , 2) , '%' )

# print("Rejected rate",int(rejected.shape)/iterations)



print("Final_Acc",final_log_acc.shape)
print("Final_rej",final_log_rej.shape)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(2,1,1)
ax.plot( rejected[:500], 'rx', label='Rejected',alpha=0.9)
ax.plot( accepted[:500], 'b.', label='Accepted',alpha=0.9)
ax.set_xlabel("Iteration")
ax.set_ylabel("v")
ax.set_title("Figure 2: MCMC sampling for v with Metropolis-Hastings. First 500 samples are shown.")
ax.grid()
ax.legend()

ax2 = fig.add_subplot(2,1,2)
to_show=final_log_acc.shape[0]
ax2.plot( final_log_rej[:to_show], 'rx', label='Rejected',alpha=0.5)
ax2.plot( final_log_acc[:to_show], 'b.', label='Accepted',alpha=0.5)
ax2.set_xlabel("Iteration")
ax2.set_ylabel("v")
ax2.set_title("Figure 3: MCMC sampling for v with Metropolis-Hastings. All samples are shown.")
ax2.grid()
ax2.legend(loc="best")

fig.tight_layout()


temp=int(0.25*all_theta.shape[0])
print("Mean_whole_theta after burn-in: ",np.mean(all_theta[temp:]))
print("SD_whole_theta after burn-in: ",np.std(all_theta[temp:]))

# consider the initial 25% of the values of v to be "burn-in", so we drop them.
show=int(0.25*final_log_acc.shape[0])

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
ax.plot(final_log_acc[show:])
ax.set_title("Figure 4: Trace for v")
ax.set_ylabel("v")
ax.set_xlabel("Iteration")





fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
# ax.hist(accepted,bins=300,label="Predicted density")
ax.set_xlabel("$\Theta$")
ax.set_ylabel("Frequency")
ax.set_title("Figure 5: posterior densities for sequence distance $\Theta$ under the JC69 model")
# ax.legend()
ax.grid("off")
ax = sns.distplot(accepted,bins=300,hist= True)

plt.show()