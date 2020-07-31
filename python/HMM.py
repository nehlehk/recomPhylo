import numpy as np
from hmmlearn import hmm
import matplotlib.pyplot as plt
import pandas as pd



df1 = pd.read_csv("/media/nehleh/295eaca0-f110-4b79-8cbe-bc99f9f61cbd/nehleh/0_Research/PhD/Data/simulationdata/recombination/clonalframe/RAxML_perSiteLLs.likelihood_GTR", sep='\s+', header=None)
ll = df1.to_numpy()
data = np.array(ll)
X = data.reshape((-1,1))
X = X[1:1000000]

mean = np.mean(X)
std = np.std(X)
# print(mean)
# print(std)


a = -2.242906769736397
astd= .001
b = -2.2394927233501036
bstd= .001



model = hmm.GaussianHMM(n_components=2, covariance_type="full" ,algorithm='viterbi' )
# model.startprob_ = np.array([0.88, 0.12])
model.startprob_ = np.array([0.98, 0.02])
model.transmat_ = np.array([[0.9999, 0.0001] , [0.0001, 0.9999]])
model.means_ = np.array([[a, astd], [b, bstd]])
model.covars_ = np.tile(np.identity(2), (2, 1, 1))


posterior = model.predict_proba(X)
# print(posterior)
print(posterior[54150:54200])
print("------------------------------------------")
print(posterior[56980:57050])
print("------------------------------------------")
print(posterior[71000:71100])

hiddenStates = model.predict(X)
# print(hiddenStates)

score = model.score(X)

fig = plt.figure(figsize=(15,8))
ax = fig.add_subplot(2,1,1)
ax.set_title("Hidden Markov Models - ClonalFrame and Recombination -- log probability of the most likely state is  " + str (score))
ax.plot(hiddenStates)
ax.set_ylabel("Clonal - NonClonal State")



ax2 = fig.add_subplot(2,1,2)
ax2.plot(posterior)
ax2.set_ylabel("posterior probability for each state")
plt.show()












# newmodel =  hmm.GaussianHMM(n_components=2 ,covariance_type="full",params = 'st').fit(X)
# newmodel.startprob_ = np.array([0.88, 0.12])
# newmodel.transmat_ = np.array([[0.9999, 0.0001] , [0.0001, 0.9999]])
# print(newmodel)
#
# print("emission-means:")
# print(newmodel.means_)
# print("emission-covar:")
# print(newmodel.covars_)
#
# print(newmodel.startprob_)
# print(newmodel.transmat_)
#
#
# posterior = newmodel.predict_proba(X)
# hiddenStates = newmodel.predict(X)
# score = newmodel.score(X)
#
# fig = plt.figure(figsize=(15,8))
# ax = fig.add_subplot(2,1,1)
# ax.set_title("Hidden Markov Models - ClonalFrame and Recombination -- log probability of the most likely state is  " + str (score))
# ax.plot(hiddenStates)
# ax.set_ylabel("Clonal - NonClonal State")
#
#
# ax2 = fig.add_subplot(2,1,2)
# ax2.plot(posterior)
# ax2.set_ylabel("posterior probability for each state")
# plt.show()

