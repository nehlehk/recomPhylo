from dendropy import Tree, DnaCharacterMatrix
import numpy
import math
import dendropy


alignment = dendropy.DnaCharacterMatrix.get(file=open("/media/nehleh/295eaca0-f110-4b79-8cbe-bc99f9f61cbd/nehleh/0_Research/PhD/Data/simulationdata/recombination/clonalframe/wholegenome.fasta"), schema="fasta")
alignment_len = alignment.sequence_size
tips = len(alignment)


column = []
for l in range(alignment_len):
    col = ""
    for t in range(tips):
        col += str(alignment[t][l])
    column.append(col)

uniqueCol = list(set(column))
# print(len(uniqueCol))
print(uniqueCol[0:10])




clonal_vector = numpy.zeros((5, 1000000))
for i in range(len(uniqueCol)):
    print(column.count(uniqueCol[i]))
    indexes = [id for id, x in enumerate(column) if x == uniqueCol[0] ]
    clonal_vector[:,indexes] = 1
    print(indexes)

print("--****************************--")
print(clonal_vector)