import numpy as np
import numpy.linalg as la
from dendropy import Tree, DnaCharacterMatrix
import myPhylo


tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_6taxa.tree'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/sample_6taxa.fasta"), schema="fasta")

tree2 = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/exampledataset/RerootTree_node12', 'newick')


pi = [0.317, 0.183 ,0.367 ,0.133]
rates = [0.000100, 0.636612 ,2.547706, 0.000100 ,2.151395]
GTR_sample = myPhylo.GTR_model(rates,pi)

column = myPhylo.get_DNA_fromAlignment(alignment)
dna = column[0]
myPhylo.set_index(tree,dna)



print("Original tree:::::::::::::::")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())

LL_normal = myPhylo.computelikelihood(tree,dna,GTR_sample)
W_LL_normal = myPhylo.wholeAlignmentLikelihood(tree,alignment,GTR_sample)



n= W_LL_normal[0]
print(n)


filter_fn = lambda n: hasattr(n, 'index') and n.index == 7
target_node = tree.find_node(filter_fn=filter_fn)
tree.reroot_at_node(target_node, update_bipartitions=False ,suppress_unifurcations = True)

LL_Reroot = myPhylo.computelikelihood(tree,dna,GTR_sample)
W_LL_rerooted = myPhylo.wholeAlignmentLikelihood(tree,alignment,GTR_sample)
myPhylo.set_index(tree,dna)

r = W_LL_rerooted[0]
print(r)

print("rooted tree:::::::::::::::")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())