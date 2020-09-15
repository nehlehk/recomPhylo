from Bio import AlignIO


path = "/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/exampledataset"

alignment = AlignIO.read(path + "wholegenome.fasta", "fasta")




a =alignment[:,1:2]
b =alignment[:,1:500]
# b = alignment[:,54176:56997]
# c = alignment[:,56997:71027]
# d = alignment[:,71027:83301]
# e = alignment[:,83301:500000]

finalAlign = b

print(finalAlign)

AlignIO.write(finalAlign , path + "sample.fasta" ,  "fasta")


# for i in range(40000,90000,5000):
#     print(i , i+5000)
#     a = alignment[:, i:i+5000]
#     AlignIO.write(a, path+ '/slices/' + str(i)+".fasta", "fasta")