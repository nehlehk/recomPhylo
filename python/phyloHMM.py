import phyloLL_emission
import pandas as pd
import numpy as np



df1 = pd.read_csv("/home/nehleh/PhyloCode/Data/RAxML_perSiteLLs_exampledatset", sep='\s+', header=None)
ll = df1.to_numpy()
data = np.array(ll)
X = data.reshape((-1,1))
X = X[1:1000000]



model = phyloLL_emission.phyloLL_HMM(n_components= 2)
