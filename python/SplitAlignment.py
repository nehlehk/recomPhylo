from Bio import AlignIO


path = "/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/20000/"

alignment = AlignIO.read(path + "wholegenome.fasta", "fasta")


# 1.   814      955      141
# 2.   1633     2247     614
# 3.   4894     5142     248
# 4.   8797     9488     691
# 5.   9488     9524     36
# 6.   9524     9530     6
# 7.   11995    12003    8
# 8.   19188    20000    812

a =alignment[:,0:814]
b =alignment[:,955:1633]
c = alignment[:,2247:4898]
d = alignment[:,5142:8797]
e = alignment[:,9530:11995]
f = alignment[:,12003:19188]
# g = alignment[:,11529:19142]

finalAlign = a + b + c + d + e + f

print(finalAlign)

AlignIO.write(finalAlign , path + "removedAllRecom.fasta" ,  "fasta")


# for i in range(40000,90000,5000):
#     print(i , i+5000)
#     a = alignment[:, i:i+5000]
#     AlignIO.write(a, path+ '/slices/' + str(i)+".fasta", "fasta")