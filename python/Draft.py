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


two_taxa = myPhylo.GTR_model(rates,pi)

c = two_taxa.get_DNA_fromAlignment(alignment)
print(c[0])
#
#
two_taxa.set_index(tree,c[0])

two_taxa.partial_likelihoods_to_target_node(tree3,c[0])







