import GTR
from dendropy import Tree, DnaCharacterMatrix
import dendropy

from GTR import GTR_model

tree = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/LL_vector/test/tree.tree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/LL_vector/test/test.fasta"), schema="fasta")

pi = [0.2184,0.2606,0.3265,0.1946]
rates = [2.0431,0.0821,0,0.067,0]

test = GTR_model(rates,pi)

res = test.computelikelihood(tree,'TTTT')

print(res[0])
print(res[1])


tmp = test.wholeAlignmentLikelihood(tree, alignment)
print(tmp[0])
print(tmp[1])

# res_exp = test.expectedLL(tree,'TTTT', 0.8, 1.2)
#
# print(res_exp[0])
# print(res_exp[1])



