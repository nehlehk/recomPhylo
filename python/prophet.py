import numpy as np
import pandas as pd


with open('/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/faketree/RAxML_perSiteLLs.likelihood_GTR', 'r') as f:
    data = f.read().split(" ")
    ll = []
    for elem in data:
        try:
            ll.append(float(elem))
        except ValueError:
            pass


signal = np.array(ll)



