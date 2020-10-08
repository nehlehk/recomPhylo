from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np
import phyloHMM
import matplotlib.pyplot as plt


def setup_indexes(tree, dna):
    """
    Set up indexes for each node of tree.
    A leaf node index is its position in the alignment.
    :param tree:
    :param dna:
    :return:
    """
    sequence_count = len(dna)

    for node in tree.postorder_node_iter():
        node.index = -1
        node.annotations.add_bound_attribute("index")

    s = sequence_count
    for node in tree.postorder_node_iter():
        # print(node.taxon)
        if not node.is_leaf():
            node.index = s
            node.label = str(node.index)
            s += 1
        else:
            for idx, name in enumerate(dna):
                # print(idx , str(name) , str(node.taxon))
                if str(name) == str(node.taxon):
                    node.index = idx
                    node.label = str(node.index)
                    break


tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/ShortDataset/RAxML_bestTree.tree'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/ShortDataset/wholegenome.fasta"), schema="fasta")


print(tree.as_ascii_plot())

pi = [0.2184,0.2606,0.3265,0.1946]
rates = [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ]
GTR_sample = myPhylo.GTR_model(rates,pi)

column = myPhylo.get_DNA_fromAlignment(alignment)
dna = column[0]
setup_indexes(tree,alignment)
tips = len(dna)


taxon = tree.taxon_namespace
nu = 0.4


# filter_fn = lambda n: hasattr(n, 'index') and n.index == 10
# target_node = tree.find_node(filter_fn=filter_fn)
# tree.reroot_at_node(target_node, update_bipartitions=False ,suppress_unifurcations = True)
# setup_indexes(tree,alignment)





def compute_logprob_phylo(X, recom_trees, model):
    n, dim = X.shape
    result = np.zeros((n, len(recom_trees)))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        children = state_tree.seed_node.child_nodes()
        for site_id, partial in enumerate(X):
            p = np.zeros(4)
            p = np.dot(model.p_matrix(children[0].edge_length), partial[0:4])
            for i in range(1, len(children)):
                p *= np.dot(model.p_matrix(children[i].edge_length), partial[i * 4:(i + 1) * 4])
            # result[site_id, tree_id] = sum(p)
            # print(p)
            site_l = np.dot(p, model.get_pi())
            result[site_id, tree_id] = np.log(site_l)
    return result






print(column[5002])

print("partial    ll old ===========================================")
persite_ll1, partial1 = myPhylo.computelikelihood(tree, column[5002],GTR_sample)
print(partial1)


# LL_root1, LL_partial1 = myPhylo.wholeAlignmentLikelihood(tree,alignment,GTR_sample)
#
# print(LL_partial1[5000])
# print(LL_root1[5000:5050])
# print(tree.as_ascii_plot(show_internal_node_labels = True))


print("partial ll new ===========================================")
myPhylo.set_index(tree,alignment)
persite_ll, partial = myPhylo.computelikelihood_new(tree, column[5002],GTR_sample)
print(partial)
# LL_root, LL_partial = myPhylo.wholeAlignmentLikelihood_new(tree,alignment,GTR_sample)
# print(LL_partial[5000])
# print(LL_root[5000:5050])
# print(LL_partial.shape)





# print("partial mixture ===========================================")
# tipdata = myPhylo.set_tips_partial(tree,alignment)
# LL_root2, LL_partial2 = myPhylo.computelikelihood_mixture(tree,alignment,tipdata,GTR_sample)
# print(LL_partial2[5000])
# print(LL_root2[5000:5050])
#
# fig = plt.figure(figsize=(15, 8))
#
# ax = fig.add_subplot(2, 1, 1)
# ax.plot(LL_root)
#
# plt.show()


