import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
import pandas as pd
from scipy.signal import blackman



df1 = pd.read_csv("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/faketree/12.txt", sep='\s+', header=None)
ll = df1.to_numpy()
data = np.array(ll)
X = data.reshape((-1,1))

w = blackman(1000)

y = fft(X*w)
print(y)

fig = plt.figure(figsize=(15,8))
plt.plot(y)
plt.show()

# print(y)
# print(y.shape)
# print(y[:,1])
# print(y[:,2])
# print(y[:,10])
# print(y[:,22500])





