from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np
import phyloLL_emission


tree = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_4taxa.tree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/sample_4taxa.fasta"), schema="fasta")


pi = [0.2184,0.2606,0.3265,0.1946]
rates = [2.0431,0.0821,0,0.067,0]


GTR_sample = myPhylo.GTR_model(rates,pi)


sitell , partial =myPhylo.wholeAlignmentLikelihood(tree,alignment, GTR_sample)


print("Before Re-rooting:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())




def reroot_tree(tree,params):
        pdm = tree.phylogenetic_distance_matrix()
        taxon = tree.taxon_namespace
        mrca = pdm.mrca(taxon[params[0]], taxon[params[1]])
        tree.reroot_at_node(mrca, update_bipartitions=False)
        return mrca


mrca = reroot_tree(tree, [0,1])
print("After Re-rooting:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())


def make_hmm_input(tree,alignment,model,params):
        reroot_tree(tree,params)
        sitell , partial =myPhylo.wholeAlignmentLikelihood(tree,alignment, model)
        children = tree.seed_node.child_nodes()
        children_count = len(children)
        x = np.zeros((alignment.sequence_size, children_count * 4))
        for id , child in enumerate(children):
                x[:,(id*4):((id+1)*4)] = partial[:, child.index, :]
        return x


# x = make_hmm_input(tree,alignment,GTR_sample,[0,2])
# print(x)
# print(x.shape)



def make_recombination_trees(tree,co_recom ,params):
    recombination_trees = []
    recombination_trees.append(tree.as_string(schema="newick"))
    tmp_tree = tree.extract_tree()
    mrca = reroot_tree(tmp_tree, params)
    print(tmp_tree.as_string(schema='newick'))
    print(tmp_tree.as_ascii_plot())
    # for node in tree.postorder_node_iter():
    #     if node.edge.length is None:
    #         node.edge.length = 0
    #     if (node.edge.length > 0) and (node.parent_node == mrca):
    #         node.edge.length = node.edge.length * co_recom
    #         recombination_trees.append(tree.as_string(schema="newick"))
    #         node.edge.length = node.edge.length * (1 / co_recom)




    # print(mrca.num_child_nodes())
    children = mrca.child_nodes()
    children[0].edge.length = children[0].edge.length + co_recom
    children[1].edge.length = children[1].edge.length + co_recom
    children[2].edge.length = children[2].edge.length - co_recom
    recombination_trees.append(tmp_tree.as_string(schema="newick"))




    mrca = reroot_tree(tree, params)
    children = mrca.child_nodes()
    children[0].edge.length = children[0].edge.length - co_recom
    children[1].edge.length = children[1].edge.length - co_recom
    children[2].edge.length = children[2].edge.length + co_recom
    recombination_trees.append(tree.as_string(schema="newick"))




    return recombination_trees


# recom_trees = make_recombination_trees(tree,0.15,[0,2])
# print(recom_trees)


# def compute_logprob_phylo(X, tree, params ,model):
#         recom_trees = make_recombination_trees(tree, 5, [0, 2])
#         n , dim = X.shape
#         result = np.zeros((n,len(recom_trees)))
#         for tree_id,item in enumerate(recom_trees):
#                 state_tree = dendropy.Tree.get(data= item, schema="newick")
#                 reroot_tree(state_tree, params)
#                 children = state_tree.seed_node.child_nodes()
#                 for site_id,partial in enumerate(X):
#                         p = np.zeros(4)
#                         p = np.dot(model.p_matrix(children[0].edge_length), partial[0:4])
#                         for i in range(1, len(children)):
#                                 p *= np.dot(model.p_matrix(children[i].edge_length), partial[i*4:(i+1)*4])
#                         result[site_id,tree_id] = sum(p)
#         return result
#
#
#
# logprob = compute_logprob_phylo(x,tree,[0,2],GTR_sample)
# print(logprob)
#
#
# model = phyloLL_emission.phyloLL_HMM(n_components= 4 ,algorithm='viterbi')

























