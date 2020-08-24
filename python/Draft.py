from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np
import phyloHMM

tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_6taxa.tree'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/sample_6taxa.fasta"), schema="fasta")



goal_tree = Tree.get_from_path('/home/nehleh/Downloads/recom_test.tree', 'newick')
# ((0:0.25,((1:0.1):0.1,2:0.2):0.05):0.05,(3:0.2,(4:0.15,5:0.15):0.05):0.1);



pi = [0.2184,0.2606,0.3265,0.1946]
rates = [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ]
GTR_sample = myPhylo.GTR_model(rates,pi)

column = myPhylo.get_DNA_fromAlignment(alignment)
myPhylo.set_index(tree,column[0])
# tmp_tree1 = tree.extract_tree()

print("Original tree")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot(plot_metric='length'))



def tree_evolver(tree ,node ,nu , position):
    recombination_trees = []
    co_recom = nu/2
    parent = node.parent_node
    print("my node:", node.index , node.edge_length)
    print("parent:" , parent.index, parent.edge_length)
    # print("sister:" , sister[0].index, sister[0].edge_length)
    if (co_recom + node.edge_length) < (node.edge_length + parent.edge_length):
        node.edge.length = node.edge.length + co_recom
        if (position == 'left') or (position == 'right'):
            sister = node.sister_nodes()
            sister[0].edge.length = sister[0].edge.length + co_recom
            parent.edge.length = parent.edge.length - co_recom
        elif (position == 'top'):
            children = node.child_nodes()
            children[0].edge.length = children[0].edge.length - co_recom
            children[1].edge.length = children[1].edge.length - co_recom
        recombination_trees.append(tree.as_string(schema="newick"))


    elif ((co_recom + node.edge_length) > (node.edge_length + parent.edge_length)) and (co_recom + node.edge_length) < tree.max_distance_from_root() :
        print("medium modification")
    elif (co_recom + node.edge_length) > tree.max_distance_from_root():
        print("final modification")

    return recombination_trees


node0 = tree.find_node_with_taxon_label(label ="0")
re0 = tree_evolver(tree, node0 , 0.1 ,'right')
print(re0)

tree = Tree.get_from_path(tree_path, 'newick')
myPhylo.set_index(tree,column[0])

node1 = tree.find_node_with_taxon_label(label ="1")
re1 = tree_evolver(tree, node1 , 0.1, 'left')
print(re1)

tree = Tree.get_from_path(tree_path, 'newick')
myPhylo.set_index(tree,column[0])

taxon = tree.taxon_namespace
pdm = tree.phylogenetic_distance_matrix()
mrca = pdm.mrca(taxon[0], taxon[1])

re2 = tree_evolver(tree, mrca , 0.1 , 'top')
print(re2)


# c = tree.nodes()
# print(len(c))
# for i,child in enumerate(c):
#     print(child.index , "::::" , child.edge_length)

# print(tree.length())
# print(tree.max_distance_from_root())
#

#
#
#
# pdm = tree.phylogenetic_distance_matrix()
# taxon = tree.taxon_namespace
# mrca = pdm.mrca(taxon[0], taxon[1])
# print("mrca:::" , mrca.index)
# p = mrca.parent_node
# print("parent:::::" , p.index)
# root = tree.seed_node
# print("root::::", root.index)
# children = mrca.child_nodes()
# print(children[0].taxon)
# print(children[1].index)

# tree.prune_subtree(children[0], update_bipartitions=True, suppress_unifurcations=True)
# mrca.remove_child(children[0])
# print("After remove:")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot(plot_metric='length'))
# print("root::::", root.index)
#
#
# dendropy.datamodel.treemodel.Node()
# taxon_1 = dendropy.Taxon(label="0")
# kid = dendropy.datamodel.treemodel.Node(edge_length= 0.03)
# root.remove_child(p)
# kid.add_child(p)
# children[0].edge_length = 0.25
# kid.add_child(children[0])
#
#
# root.add_child(kid)
# # children[0].edge_length = 0.3
# tree.ladderize(ascending=False)
#
# print("After add:")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot(plot_metric='length'))



# mrca = myPhylo.reroot_tree(tree,[0,1])
# adjacent_childs = tree.seed_node.child_nodes()
# childA = adjacent_childs[0]
# childB = adjacent_childs[1]
# childC = adjacent_childs[2]
# print("After Re-rooting:")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot(plot_metric='length'))
#
#
# print(root.index)
# tree.reroot_at_node(root, update_bipartitions=False)
# print("Original tree?:")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot(plot_metric='length'))



# mrca.parent_node = tree.seed_node
# print(mrca.edge_length)
# mrca.edge_length = 0.5
# nane = mrca.parent_node
# print(nane.index)
#
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot(plot_metric='length'))


# tree.reseed_at(mrca,update_bipartitions=True, collapse_unrooted_basal_bifurcation=True , suppress_unifurcations=True)
# mrca = myPhylo.reroot_tree(tree,[0,1])
# print("After Re-seed:")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot(plot_metric='length'))
#
#
# childs = tree.seed_node.child_nodes()
# print(tree.seed_node.index)
# print(childs[0].index)
# print(childs[1].index)
# print(childs[2].index)
#
# nane = childs[2].parent_node
#
# print("parent::::" , nane.index)
#
#
# tree.prune_subtree(childs[0] ,update_bipartitions=True, suppress_unifurcations=False)
#
#
#
#
# print("final tree:")
# print(goal_tree.as_string(schema='newick'))
# print(goal_tree.as_ascii_plot(plot_metric='length'))


# mrca.remove_child( childs[2] ,suppress_unifurcations=False)
# print(mrca.parent_node())
# mrca2 = pdm.mrca(taxon[1], taxon[2])
#
# tree.seed_node.new_child()

# print("After prune:")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot())


# tmp_tree1 = tree.extract_tree()
# mrca = myPhylo.reroot_tree(tree, [0,1])
# print(mrca.index)
# children = mrca.child_nodes()
# for i in range(mrca.num_child_nodes()):
#     print(children[i].index)
# tree.prune_taxa('0' ,update_bipartitions=True ,is_apply_filter_to_leaf_nodes=True,)
# mrca.remove_child(children[0] ,suppress_unifurcations=True)

# mrca.insert_new_child()


# print("After prune:")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot())
#
# print(mrca.parent_node)


# x =  myPhylo.make_hmm_input(tree,alignment,GTR_sample,[0,2])
# print(x)



# recom_trees = myPhylo.make_recombination_trees(tree,0.05,[0,1])
# print(recom_trees)































