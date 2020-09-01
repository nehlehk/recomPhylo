import numpy as np
import numpy.linalg as la
from dendropy import Tree, DnaCharacterMatrix
import myPhylo


tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/exampledataset/exampledataset_RAxML_bestTree'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/exampledataset/wholegenome.fasta"), schema="fasta")


pi = [0.2184,0.2606,0.3265,0.1946]
rates = [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ]
GTR_sample = myPhylo.GTR_model(rates,pi)

column = myPhylo.get_DNA_fromAlignment(alignment)
myPhylo.set_index(tree,column[0])
dna = column[0]


print("Original tree:::::::::::::::")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())

# LL_normal = myPhylo.computelikelihood(tree,dna,GTR_sample)
# W_LL_normal = myPhylo.wholeAlignmentLikelihood(tree,alignment,GTR_sample)
#
#
# print("LL_normal:" , LL_normal[0])
# n= W_LL_normal[0]
# print(n[0:25])

filter_fn = lambda n: hasattr(n, 'index') and n.index == 14
target_node = tree.find_node(filter_fn=filter_fn)
tree.reroot_at_node(target_node, update_bipartitions=False ,suppress_unifurcations = True)

# LL_Reroot = myPhylo.computelikelihood(tree,dna,GTR_sample)
# W_LL_rerooted = myPhylo.wholeAlignmentLikelihood(tree,alignment,GTR_sample)
#
# print("LL_Reroot:" , LL_Reroot[0])
# r = W_LL_rerooted[0]
# print(r[0:25])

print("Re_root tree:::::::::::::::")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())