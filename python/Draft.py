from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np


tree = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_2taxa.tree', 'newick')
tree2 = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_2taxa2.tree', 'newick')
tree3 = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_2taxa3.tree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/sample_2taxa.fasta"), schema="fasta")




pi = [0.2184,0.2606,0.3265,0.1946]
rates = [2.0431,0.0821,0,0.067,0]


# phylo_m.traverse()


two_taxa = myPhylo.GTR_model(rates,pi)
#
c = two_taxa.get_DNA_fromAlignment(alignment)
print(c[0])
#
#
two_taxa.set_index(tree,c[0])



# peeling = two_taxa.get_peeling_order(tree3)
# print(peeling)

#
#
a1 = two_taxa.computelikelihood(tree, c[2])
# print(a1[0])
# print(a1[1])
#
# a2 = two_taxa.computelikelihood(tree2,c[2])
# # print(a2[0])
# #
# a3 = two_taxa.computelikelihood(tree3,c[2])
# # print(a3[0])
#
#
#
# tmp = two_taxa.wholeAlignmentLikelihood(tree, alignment)
# print(tmp[0])
#
# tmp2 = two_taxa.wholeAlignmentLikelihood(tree2, alignment)
# print(tmp2[0])
#
# tmp3 = two_taxa.wholeAlignmentLikelihood(tree3, alignment)
# print(tmp3[0])
#




