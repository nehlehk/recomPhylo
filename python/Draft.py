from dendropy import Tree, DnaCharacterMatrix
import dendropy
from GTR import GTR_model

tree = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/exampledataset/tree.tree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/exampledataset/wholegenome.fasta"), schema="fasta")

pi = [0.2184,0.2606,0.3265,0.1946]
rates = [2.0431,0.0821,0,0.067,0]

test = GTR_model(rates,pi)

c = test.get_DNA_fromAlignment(alignment)
print(c[100])

# res = test.computelikelihood(tree, c[0])
# print(res[0])
# print(res[1])


tmp = test.wholeAlignmentLikelihood(tree, alignment)
print(tmp[0].shape)
print(tmp[1].shape)

# res_exp = test.wholeAlignmentExpLL(tree,alignment, 0.8, 1.2)
# print(res_exp[0])
# print(res_exp[1])



