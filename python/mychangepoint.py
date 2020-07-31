# import libraries
import numpy as np
import matplotlib.pyplot as plt


blocks = np.loadtxt(open("/home/nehleh/0_Research/PhD/Data/xmfa/xmfa-blocks.csv", "rb"), delimiter=",", skiprows=1)

myblock = blocks[:,1]




# creation of data
with open('/home/nehleh/0_Research/PhD/Data/Ecoli-10_Aligned/likelihood_persite_JC', 'r') as f:
    data = f.read().split(" ")
    ll = []
    for elem in data:
        try:
            ll.append(float(elem))
        except ValueError:
            pass


signal = np.array(ll)

alignlen = len(signal)

mean = np.mean(signal)
std = np.std(signal)

maxdata= max(signal)
mindata= min(signal)

# print(mean)
# print(maxdata)
# print(mindata)

x = np.arange(0, alignlen, 1)

# less = np.ma.masked_where(signal > mean - 2*std , signal)
# middle = np.ma.masked_where((signal <  mean + 2*std ) & (signal >   mean - 2*std ), signal)
# more = np.ma.masked_where(signal < 5* mean , signal)


plt.figure(figsize=(25,10))
# plt.plot( x, more , x,less , x , middle)
# plt.plot(x, signal , color = "g")
plt.plot( x , signal)
for i in (myblock):
    plt.axvline(i, color='r', ls='-')
plt.show()
