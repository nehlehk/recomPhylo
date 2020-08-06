from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np


tree = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_4taxa.tree', 'newick')
tree2 = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_2taxa2.tree', 'newick')
tree3 = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_2taxa3.tree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/sample_2taxa.fasta"), schema="fasta")


pi = [0.2184,0.2606,0.3265,0.1946]
rates = [2.0431,0.0821,0,0.067,0]


two_taxa = myPhylo.GTR_model(rates,pi)

# c = two_taxa.get_DNA_fromAlignment(alignment)
# print(c[0])
#
#
# two_taxa.set_index(tree,c[0])

# node_partial=two_taxa.partial_likelihoods_to_target_node(tree,'AAAA')
#
# print(node_partial)


def make_recombination_tree(tree ,co_recom):
        recombination_trees = []
        for node in tree.postorder_node_iter():
            if node.edge.length is None:
                node.edge.length = 0
            # print(node.edge.length)
            if node.edge.length > 0:
                node.edge.length = node.edge.length * co_recom
                recombination_trees.append(tree.as_string(schema="newick"))
                node.edge.length = node.edge.length * 1/co_recom

        return recombination_trees




recom = myPhylo.make_recombination_trees(tree3 ,8)
print(recom[1])






