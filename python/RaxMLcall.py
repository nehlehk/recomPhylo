import os


# bootstrapp = " -m GTRCAT -p 12345 -b 12345  -V -# 100 -s /home/nehleh/0_Research/PhD/Data/LL_vector/JC69_100.fasta  -n T14"
# bipartition = " -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T13 -z RAxML_bootstrap.T14 -n T15"

path = "/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/exampledataset/"

help = "/home/nehleh/Documents/0_Research/Software/standard-RAxML-master/raxmlHPC -h"

raxml = "/home/nehleh/Documents/0_Research/Software/standard-RAxML-master/raxmlHPC"

param = " -m GTRGAMMA   -p 12345 -s " + path+ "wholegenome.fasta -N 10 -n wholegenometree"

param2 = " -m GTRGAMMA  --JC69  -p 12345 -s " + path+ "wholegenome.fasta -n likelihood_JC "

paramlikelihood = " -m GTRGAMMA   -p 12345 -s " + path+ "wholegenome.fasta -n likelihood_GTR_Reroot14"

likelihood = " -f g -z " + path+ "RerootTree_node14"

cmd = raxml + paramlikelihood + likelihood
# cmd = raxml + param
print(cmd)
os.system(cmd)