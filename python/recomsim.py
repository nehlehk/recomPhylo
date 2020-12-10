import os

seq_count = '10'
seqnum = '5000'
# =====================================  fastsimbactree ====================================

# fastsimbac = "/home/nehleh/Documents/anaconda3/pkgs/fastsimbac-1.0.1_bd3ad13d8f79-h6dcb523_0/bin/fastSimBac  "
#
# # paramfast = seq_count + "  " + seqnum + " -T -t .001 -r 0 0  -x 0.0005 100  −I 20 100 0.0 −ej 20.0 1 2 > /home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/test/fastsimbactrees.txt"
# paramfast = seq_count + "  " + seqnum + " -T -t .001 -r 0 0 -x 0.0001 1000 -I 2 0 6 -ej 20 1 2 -d > /home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/test/fastsimbactrees.txt"
#
# cmdfastsimbac = fastsimbac + paramfast
# print(cmdfastsimbac)
# os.system(cmdfastsimbac)

# =====================================  prepare tree for seq-gen ====================================
# f = open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/test/fastsimbactrees.txt", "r")
treefile= "/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/test/tree.tree"
partition = 11
# for line in f:
#     if  line.find('[') == 0:
#         partition += 1
#         fi = open(treefile, "a")
#         fi.write(line)
#         fi.close()

# =====================================  seqgen ====================================
model = 'GTR'

frequencies = '0.2184,0.2606,0.3265,0.1946'
rates = '0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1'
# frequencies = '0.25,0.25,0.25,0.25'
# rates = '1,1,1,1,1,1'
outfile = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/test/wholegenome.fasta'
seqgen = '/home/nehleh/Documents/0_Research/Software/Seq-Gen-1.3.4/source/seq-gen'
cmd = seqgen + '  -m'+model+ '  -l'+ seqnum + '  -f'+frequencies  + '  -p' + str(partition)  \
      +  '  -s0.2' + '  -r'+rates +'  -of'+ '  <'+treefile+'>  ' + outfile
print(cmd)
os.system(cmd)


