from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np


# ==============================================   input  ==============================================================
tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_6taxa.tree'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/sample_6taxa.fasta"), schema="fasta")


pi = [0.2184,0.2606,0.3265,0.1946]
rates = [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ]
GTR_sample = myPhylo.GTR_model(rates,pi)

column = myPhylo.get_DNA_fromAlignment(alignment)
myPhylo.set_index(tree,column[0])


taxon = tree.taxon_namespace
# ============================================  methods ================================================================
def find_kids_index(node):
    kids = []
    for id, child in enumerate(node.child_node_iter()):
        kids.append(child.index)
    return kids
# ----------------------------------------------------------------------------------------------------------------------
def find_sibling_index(node):
    s = []
    sibling = node.sibling_nodes()
    for i in range(len(sibling)):
        s.append(sibling[i].index)
    return s
# ----------------------------------------------------------------------------------------------------------------------
def find_node_position(node,target_kids):
    sister = find_sibling_index(node)
    if sister[0] in target_kids:
        return "descendant"
    else:
        return "ancestor"
# ----------------------------------------------------------------------------------------------------------------------


# ==============================================   input  ==============================================================
print("Original tree:::::::::::::::")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())




# mytree = []
#
# for id_tree, target_node in enumerate(tree.postorder_internal_node_iter(exclude_seed_node=True)):
#     print(target_node.index)
#     # print(id_tree)
#     mytree.append(Tree.get_from_path(tree_path, 'newick'))
#     myPhylo.set_index(mytree[id_tree], column[0])
#     # ----------- Step 1 : Make input for hmm ------------------------------------------------------
#     # --------------  Stetp 1.1 : We need to re-root the tree based on the target node where the target node is each internal node of the tree.
#     mytree[id_tree].reroot_at_node(target_node, update_bipartitions=False ,suppress_unifurcations = True)
#
#     # --------------  Step 1.2: Calculate X based on this re-rooted tree
#     X = myPhylo.make_hmm_input(mytree[id_tree], alignment, GTR_sample)
#     # print(X)
#
#     # --------------  Step 1.3: Reset the rerooted to the original tree
#     tree = Tree.get_from_path(tree_path, 'newick')
#     myPhylo.set_index(tree, column[0])
#     mytree[id_tree] = Tree(tree)
#
#     # ----------- Step 2: make 3 recombination trees -----------------------------------------------()
#     recombination_trees = []
#     temptree = {}
#     target_kids = find_kids_index(target_node)
#     for id, child in enumerate(target_kids):
#         print(child)
#         temptree["tree{}".format(id)] = Tree.get_from_path(tree_path, 'newick')
#         myPhylo.set_index(temptree["tree{}".format(id)], column[0])
#         filter_fn = lambda n: hasattr(n, 'index') and n.index == child
#         recombined_node = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
#         position = find_node_position(recombined_node,target_kids)
#         print(position)
#         tmp = myPhylo.tree_evolver(temptree["tree{}".format(id)],recombined_node,.24 , position)
#         print(tmp)

    # for id, child in enumerate(target_node.child_node_iter()):
    #     print(child.index)
    #     temptree["tree{}".format(id)] = Tree.get_from_path(tree_path, 'newick')
    #     myPhylo.set_index(temptree["tree{}".format(id)], column[0])
    #     i = child.index
    #     filter_fn = lambda n: hasattr(n, 'index') and n.index == i
    #     recombined_node = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
    #     tmp = myPhylo.tree_evolver(temptree["tree{}".format(id)],recombined_node,.24)
    #     print(tmp)

    # ----------- Step 3: Call phyloHMM ----------------------------------------------------------





# ======================================================================================================================
filter_fn = lambda n: hasattr(n, 'index') and n.index == 8
target_node = tree.find_node(filter_fn=filter_fn)
tree.reroot_at_node(target_node, update_bipartitions=False ,suppress_unifurcations = True)


for id, child in enumerate(target_node.child_node_iter()):
    print(child.index)
    print(child.edge_length)

# temptree = {}
# target_kids = find_kids_index(target_node)
# for id, child in enumerate(target_kids):
#     print(child)
#     temptree["tree{}".format(id)] = Tree.get_from_path(tree_path, 'newick')
#     myPhylo.set_index(temptree["tree{}".format(id)], column[0])
#     filter_fn = lambda n: hasattr(n, 'index') and n.index == child
#     recombined_node = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
#     position = find_node_position(recombined_node, target_kids)
#     print(position)
#     tmp = myPhylo.tree_evolver(temptree["tree{}".format(id)], recombined_node, .24, position)
#     print(tmp)
