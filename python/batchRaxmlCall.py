import os


fasta = os.listdir("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/500000-1/slices")

path = "/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/500000-1/slices/"

raxml = "/home/nehleh/0_Research/Software/standard-RAxML-master/raxmlHPC"

# param = " -m GTRGAMMA   -p 12345 -s " + path+ fasta[i] + " -N 10 -n "+fasta[i]

for i in range(len(fasta)):
    cmd = raxml + " -m GTRGAMMA   -p 12345 -s " + path+ fasta[i] + " -N 10 -n tree" + fasta[i].strip('.fasta')
    print(cmd)
    os.system(cmd)

# cmd = raxml + param
# print(cmd)
# os.system(cmd)