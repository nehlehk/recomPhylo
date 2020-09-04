from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np
import matplotlib.pyplot as plt


# ==============================================   input  ==============================================================
tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_6taxa.tree'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/sample_6taxa.fasta"), schema="fasta")


pi = [0.317, 0.183 ,0.367 ,0.133]
rates = [0.000100, 0.636612 ,2.547706, 0.000100 ,2.151395]

GTR_sample = myPhylo.GTR_model(rates,pi)

column = myPhylo.get_DNA_fromAlignment(alignment)
dna = column[6]
myPhylo.set_index(tree,dna)

print("Original tree:::::::::::::::")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())

# persite_ll, partial = myPhylo.computelikelihood(tree,dna,GTR_sample)

# print(persite_ll)
# print(partial)


myPhylo.computelikelihood_2(tree,dna,GTR_sample)


filter_fn = lambda n: hasattr(n, 'index') and n.index == 7
target_node = tree.find_node(filter_fn=filter_fn)
tree.reroot_at_node(target_node, update_bipartitions=False ,suppress_unifurcations = True)
myPhylo.set_index(tree,dna)


# R_persite_ll, R_partial  = myPhylo.computelikelihood(tree,dna,GTR_sample)
myPhylo.computelikelihood_2(tree,dna,GTR_sample)

#
# print(R_persite_ll)
# print(R_partial)


# print("Re_root tree:::::::::::::::")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot())



