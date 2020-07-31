import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df1 = pd.read_csv("https://raw.githubusercontent.com/nehlehk/FirstPhylo/master/Data/RAxML_perSiteLLs_exampledatset", sep='\s+', header=None)
print(df1)
ll = df1.to_numpy()
print(ll)
data = np.array(ll)
X = data.reshape((-1,1))
print(X)

fig = plt.figure(figsize=(8, 14))
sns.set_style('whitegrid')
ax = fig.add_subplot(4,1,1)
ax.set_title("Value for X[0:22000] and X[22000:25000]")
sns.distplot(X[0:22000] , label= 'Clonal'  )
sns.distplot(X[22000:25000] , label= 'recombination' )
ax.legend()


ax2 = fig.add_subplot(4,1,2)
ax2.set_title("Value for X[25000:47000] and X[47000:50000]")
sns.distplot(X[25000:47000]  , label= 'Clonal' )
sns.distplot(X[47000:50000] , label= 'recombination' )
ax2.legend()

ax3 = fig.add_subplot(4,1,3)
ax3.set_title("Value for X[50000:72000] and X[72000:75000]")
sns.distplot(X[50000:72000]  , label= 'Clonal' )
sns.distplot(X[72000:75000] , label= 'recombination' )
ax3.legend()

ax4 = fig.add_subplot(4,1,4)
ax4.set_title("Value for X[75000:97000] and X[97000:100000]")
sns.distplot(X[75000:97000]  , label= 'Clonal' )
sns.distplot(X[97000:100000] , label= 'recombination' )
ax4.legend()


plt.show()