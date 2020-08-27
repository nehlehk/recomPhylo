from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo

# ==============================================   input  ==============================================================
tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/exampledataset/tree.tree'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/exampledataset/wholegenome.fasta"), schema="fasta")


pi = [0.2184,0.2606,0.3265,0.1946]
rates = [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ]
GTR_sample = myPhylo.GTR_model(rates,pi)

column = myPhylo.get_DNA_fromAlignment(alignment)
myPhylo.set_index(tree,column[0])

# ============================================  methods ================================================================
def reset_tree():
    tree = Tree.get_from_path(tree_path, 'newick')
    myPhylo.set_index(tree, column[0])
# *************************************************************


# ==============================================   input  ==============================================================
# print("Original tree")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot(plot_metric='length'))


internl_nodes = tree.postorder_internal_node_iter(exclude_seed_node=True)
target_nodes = []
for node in internl_nodes:
    target_nodes.append(node)

for i in range(len(target_nodes)):
    print(target_nodes[i].index)



# for i in range(len(target_nodes)):
#     print("Original tree:::::::::::::::")
#     print(tree.as_string(schema='newick'))
#     print(tree.as_ascii_plot())
#     print("target node is ", target_nodes[i].index)
#     tree.reroot_at_node(target_nodes[i], update_bipartitions = True)
#     print("Reroot tree:::::::::::::::")
#     print(tree.as_string(schema='newick'))
#     print(tree.as_ascii_plot())
#     # print("seed node after reroot:", tree.seed_node.index)
#     # print(tree.as_string(schema='newick'))
#     # children = tree.seed_node.child_nodes()
#     # for id, child in enumerate(children):
#     #     print(child.index)
#     X = myPhylo.make_hmm_input(tree,alignment,GTR_sample,target_nodes[i])
#     print("After hmm tree:::::::::::::::")
#     print(tree.as_string(schema='newick'))
#     print(tree.as_ascii_plot())
#     # print(X)
#     tree = Tree.get_from_path(tree_path, 'newick')
#     myPhylo.set_index(tree, column[0])
#     print("After reset:::::::::::::::")
#     print(tree.as_string(schema='newick'))
#     print(tree.as_ascii_plot())
#     print("********************************************************************************************************")



