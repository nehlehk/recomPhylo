from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np


tree = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_2taxa.tree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/sample_2taxa.fasta"), schema="fasta")


pi = [0.2184,0.2606,0.3265,0.1946]
rates = [2.0431,0.0821,0,0.067,0]


two_taxa = myPhylo.GTR_model(rates,pi)

align = two_taxa.get_DNA_fromAlignment(alignment)

# two_taxa.set_index(tree,c[0])
recom= myPhylo.make_recombination_trees(tree ,5)
print(recom)
# print(len(recom))

x = np.zeros((len(align) , 2 * len(recom)))
# print(x.shape)
for idx,col  in enumerate(align):
    for rid,rtree in enumerate(recom):
        recom_tree = Tree.get_from_string(rtree,'newick')
        node_partial = two_taxa.partial_likelihoods_to_target_node(recom_tree,col)
        x[idx][rid] = node_partial[0]
        x[idx][rid+len(recom)] = node_partial[1]



print(x)














