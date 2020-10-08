import os


# bootstrapp = " -m GTRCAT -p 12345 -b 12345  -V -# 100 -s /home/nehleh/0_Research/PhD/Data/LL_vector/JC69_100.fasta  -n T14"
# bipartition = " -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T13 -z RAxML_bootstrap.T14 -n T15"

path = "/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/oneMilion/"

help = "/home/nehleh/Documents/0_Research/Software/standard-RAxML-master/raxmlHPC -h"

raxml = "/home/nehleh/Documents/0_Research/Software/standard-RAxML-master/raxmlHPC"

param = " -m GTRGAMMA   -p 12345 -s " + path+ "wholegenome.fasta -N 10 -n tree"

param2 = " -m GTRGAMMA  --JC69  -p 12345 -s " + path+ "wholegenome.fasta -n likelihood_JC "

paramlikelihood = " -m GTRGAMMA   -p 12345 -s " + path+ "sample_6taxa.fasta -n 6taxa_edgeOne"

likelihood = " -f g -z " + path+ "tree_6taxa_edgeone.tree"

# cmd = raxml + paramlikelihood + likelihood
cmd = raxml + param
print(cmd)
os.system(cmd)