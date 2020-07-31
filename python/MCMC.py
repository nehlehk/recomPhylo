from builtins import print

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats




population = np.random.normal(10, 3, 30000)
observation = population[np.random.randint(0, 30000, 1000)]

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1)
ax.hist(observation, bins=35)
ax.set_xlabel("value")
ax.set_ylabel("Frequency")
ax.set_title("Figure 1: Distribution of 1000 observations sampled from a population of 30,000 with mu=10, sigma=3")
mu_obs = observation.mean()
print(mu_obs)



#=======================================================================================================================
#The tranistion model defines how to move from sigma_current to sigma_new
def transition_model(x):
   # return  (x[0] ,np.random.normal(x[1], 0.5, (1,)))
   return (np.random.normal(x[0], 0.5, (1,)), np.random.normal(x[1], 0.5, (1,)))
#=======================================================================================================================
def prior(x):
    if x[1] <= 0:
        return 0
    else:
        return 1
#=======================================================================================================================
#Computes the likelihood of the data given a sigma (new or current) according to equation (2)
def manual_log_like_normal(x,data):
    # x[0]=mu, x[1]=sigma (new or current)
    # data = the observation
    mll_n = np.sum(-np.log(x[1] * np.sqrt(2* np.pi) )-((data-x[0])**2) / (2*x[1]**2))
    return mll_n
# =======================================================================================================================
#Same as manual_log_like_normal(x,data), but using scipy implementation. It's pretty slow.
def log_lik_normal(x,data):

    # x[0]=mu, x[1]=sigma (new or current)
    # data = the observation
    lln = np.sum(np.log(scipy.stats.norm(x[0],x[1]).pdf(data)))
    return lln
# =======================================================================================================================
#Defines whether to accept or reject the new sample
def acceptance(x, x_new):
    if x_new > x:
        return True
    else:
        accept = np.random.uniform(0,1)
        # Since we did a log likelihood, we need to exponentiate in order to compare to the random number
        # less likely x_new are less likely to be accepted
        return  (accept < (np.exp(x_new-x)))
# =======================================================================================================================
def metropolis_hastings(likelihood_computer,prior, transition_model, param_init,iterations,data,acceptance_rule):
    # likelihood_computer(x,data): returns the likelihood that these parameters generated the data
    # transition_model(x): a function that draws a sample from a symmetric distribution and returns it
    # param_init: a starting sample
    # iterations: number of accepted to generated
    # data: the data that we wish to model
    # acceptance_rule(x,x_new): decides whether to accept or reject the new sample


    x = param_init
    accepted = []
    rejected = []
    for i in range(iterations):
        x_new = transition_model(x)
        x_lik = likelihood_computer(x,data)
        x_new_lik = likelihood_computer(x_new,data)
        if (acceptance(x_lik+np.log(prior(x)),x_new_lik+np.log(prior(x_new)))):
            x = x_new
            accepted.append(x_new)
        else:
            rejected.append(x_new)


    return np.array(accepted) , np.array(rejected)
# =======================================================================================================================
accepted, rejected =  metropolis_hastings(manual_log_like_normal,prior,transition_model,[8,0.1],50000,observation,acceptance)

print("Accepeted",accepted.shape)

print("Rejected",rejected.shape)

# ================================plot $\sigma$ =========================
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(2,1,1)
ax.plot( rejected[0:50,1], 'rx', label='Rejected',alpha=0.9)
ax.plot( accepted[0:50,1], 'b.', label='Accepted',alpha=0.9)
ax.set_xlabel("Iteration")
ax.set_ylabel("$\sigma$")
ax.set_title("Figure 2: MCMC sampling for $\sigma$ with Metropolis-Hastings. First 50 samples are shown.")
ax.grid()
ax.legend()


ax2 = fig.add_subplot(2,1,2)
to_show=accepted.shape[0]
ax2.plot( rejected[:to_show,1], 'rx', label='Rejected',alpha=0.5)
ax2.plot( accepted[:to_show,1], 'b.', label='Accepted',alpha=0.5)
ax2.set_xlabel("Iteration")
ax2.set_ylabel("$\sigma$")
ax2.set_title("Figure 3: MCMC sampling for $\sigma$ with Metropolis-Hastings. All samples are shown.")
ax2.grid()
ax2.legend()

fig.tight_layout()
# ================================ plot $\mu$ =========================
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(2,1,1)
ax.plot( rejected[0:50,0], 'rx', label='Rejected',alpha=0.9)
ax.plot( accepted[0:50,0], 'g.', label='Accepted',alpha=0.9)
ax.set_xlabel("Iteration")
ax.set_ylabel("$\mu$")
ax.set_title("Figure 4: MCMC sampling for $\mu$ with Metropolis-Hastings. First 50 samples are shown.")
ax.grid()
ax.legend()


ax2 = fig.add_subplot(2,1,2)
to_show=accepted.shape[0]
ax2.plot( rejected[:to_show,0], 'rx', label='Rejected',alpha=0.5)
ax2.plot( accepted[:to_show,0], 'g.', label='Accepted',alpha=0.5)
ax2.set_xlabel("Iteration")
ax2.set_ylabel("$\mu$")
ax2.set_title("Figure 5: MCMC sampling for $\mu$ with Metropolis-Hastings. All samples are shown.")
ax2.grid()
ax2.legend()

fig.tight_layout()



# ========================================================================
# We consider the initial 25% of the values of $\sigma$ to be "burn-in", so we drop them.
show=int(0.25*accepted.shape[0])

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
ax.plot(accepted[show:,1])
ax.set_title("Figure 6: Trace for $\sigma$")
ax.set_ylabel("$\sigma$")
ax.set_xlabel("Iteration")
# ========================================================================
# We consider the initial 25% of the values of $\mu$ to be "burn-in", so we drop them.
show=int(0.25*accepted.shape[0])

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
ax.plot(accepted[show:,0])
ax.set_title("Figure 7: Trace for $\mu$")
ax.set_ylabel("$\mu$")
ax.set_xlabel("Iteration")
# ========================================================================
# ================================ comparison ====================
mu=accepted[show:,0].mean()
sigma=accepted[show:,1].mean()
print(mu, sigma)
model = lambda t,mu,sigma:np.random.normal(mu,sigma,t)
observation_gen=model(population.shape[0],mu,sigma)
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
ax.hist( observation_gen,bins=70 ,label="Predicted distribution of 30,000 individuals")
ax.hist( population,bins=70 ,alpha=0.5, label="Original values of the 30,000 individuals")
ax.set_xlabel("Mean")
ax.set_ylabel("Frequency")
ax.set_title("Figure 5: Posterior distribution of predicitons")
ax.legend()
ax.grid("off")


plt.show()