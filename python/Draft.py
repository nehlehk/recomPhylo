from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np


tree = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_6taxa.tree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/sample_6taxa.fasta"), schema="fasta")


pi = [0.2184,0.2606,0.3265,0.1946]
rates = [2.0431,0.0821,0,0.067,0]


GTR_sample = myPhylo.GTR_model(rates,pi)


print("Before Re-rooting:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())

def make_hmm_input(tree,alignment,model,params):
        pdm = tree.phylogenetic_distance_matrix()
        taxon = tree.taxon_namespace
        mrca = pdm.mrca(taxon[params[0]], taxon[params[1]])
        tree.reroot_at_node(mrca, update_bipartitions=False)
        sitell , partial =myPhylo.wholeAlignmentLikelihood(tree,alignment, model)
        children = tree.seed_node.child_nodes()
        children_count = len(children)
        x = np.zeros((alignment.sequence_size, children_count * 4))
        for id , child in enumerate(children):
                x[:,(id*4):((id+1)*4)] = partial[:, child.index, :]
        return x


# x = make_hmm_input(tree,alignment,GTR_sample,[0,7])
# print(x)
# print(x.shape)



def make_recombination_trees(tree,alignment,co_recom ,params):
    myPhylo.set_index(tree,alignment)
    recombination_trees = []
    recombination_trees.append(tree.as_string(schema="newick"))
    pdm = tree.phylogenetic_distance_matrix()
    taxon = tree.taxon_namespace
    mrca = pdm.mrca(taxon[params[0]], taxon[params[1]])
    tree.reroot_at_node(mrca, update_bipartitions=False)
    for node in tree.postorder_node_iter():
        if node.edge.length is None:
            node.edge.length = 0
        # print(node.edge.length)
        if (node.edge.length > 0) :
            node.edge.length = node.edge.length * co_recom
            recombination_trees.append(tree.as_string(schema="newick"))
            node.edge.length = node.edge.length * (1 / co_recom)

    return recombination_trees


recom_trees = make_recombination_trees(tree,alignment,5,[0,2])
print(recom_trees)









# align = myPhylo.get_DNA_fromAlignment(alignment)
#
# res =myPhylo.wholeAlignmentLikelihood(tree,alignment, GTR_sample)
# print(res[0])
# print(res[1])
#
#
# myPhylo.set_index(tree,align[0])
#
#
# print("seed before reroot: ",tree.seed_node.index)
#
#
# print("Before:")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot())
#
# pdm = tree.phylogenetic_distance_matrix()
# taxon = tree.taxon_namespace
# mrca = pdm.mrca(taxon[0] , taxon[2])

# print(taxon)
# print("mrca: ",mrca.index)
# print(mrca.edge_length)
# print(mrca.description())

# tree.reroot_at_node(mrca, update_bipartitions=False)
#
# print("After:")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot())
#
# print("seed after reroot: ",tree.seed_node.index)
# sitell , partial =  myPhylo.computelikelihood(tree,align[8],GTR_sample)
# print(sitell)
# print(partial)
# children = tree.seed_node.child_nodes()
# for child in children:
#         # print(child.index)
#         print(partial[child.index])























