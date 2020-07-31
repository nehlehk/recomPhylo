from statsmodels.tsa.ar_model import AutoReg
from random import random
import matplotlib.pyplot as plt
from statsmodels.tsa.holtwinters import SimpleExpSmoothing
import pandas as pd
import numpy as np


df1 = pd.read_csv("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/exampledataset/12.txt", sep='\s+', header=None)
ll = df1.to_numpy()
data = np.array(ll)
X = data.reshape((-1,1))

# # fit model
model = SimpleExpSmoothing(X)
model_fit = model.fit()
# # make prediction
yhat = model_fit.predict(len(X), len(X))
print(yhat)
plt.plot(X)
plt.axvline(yhat, color='r', ls='-.')
plt.show()


# contrived dataset
# data = [x + random() for x in range(1, 100)]
# print(data)
# plt.plot(data)
#
# # fit model
# model = AutoReg(data, lags=1)
# model_fit = model.fit()
# # make prediction
# yhat = model_fit.predict(len(data), len(data))
# print(yhat)
# plt.axvline(yhat, color='r', ls='-.')
# plt.show()