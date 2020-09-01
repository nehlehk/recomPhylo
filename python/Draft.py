from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np
import phyloHMM

tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_6taxa.tree'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/sample_6taxa.fasta"), schema="fasta")




pi = [0.2184,0.2606,0.3265,0.1946]
rates = [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ]
GTR_sample = myPhylo.GTR_model(rates,pi)

column = myPhylo.get_DNA_fromAlignment(alignment)
myPhylo.set_index(tree,column[0])

print("Original tree")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot(plot_metric='length'))



def tree_evolver2(tree ,node ,nu):
    recombination_trees = []
    co_recom = nu/2
    parent = node.parent_node
    grandparent = parent.parent_node
    br_length = node.edge_length
    print("My node is:" , node.index , node.edge_length)
    print("parent::" , parent.index)
    print("grandparent::", grandparent.index)
    print(grandparent.distance_from_tip())


    # recombination_trees = []
    # recombination_tree = ''
    # co_recom = nu/2
    # parent = node.parent_node
    # grandparent = parent.parent_node
    # br_length = node.edge_length
    # print("My node is:" , node.index , node.edge_length)
    # print("parent::" , parent.index)
    # print("grandparent::", grandparent.index)
    # print(grandparent.distance_from_tip())
    # print(tree.as_ascii_plot())
    # print(tree.as_string(schema='newick'))

    # topology does not change in this case:
    if (co_recom + br_length) < (grandparent.distance_from_tip()):
        print(" *********** Stage one ***********")
        node.edge_length = br_length + co_recom
        sister = node.sister_nodes()
        sister[0].edge_length = sister[0].edge_length + co_recom
        parent.edge_length = parent.edge_length - co_recom
        print(parent.edge_length)
        recombination_trees.append(tree.as_string(schema="newick"))
        # recombination_tree = tree.as_string(schema="newick")

    # changing in topology to make recombination tree:
    elif ((co_recom + br_length) > grandparent.distance_from_tip())  and ((co_recom + br_length) < tree.max_distance_from_root()):
        print(" *********** Stage Two ***********")
        ancestor = []
        recom_length = co_recom + node.edge_length
        for id,tmp_node in enumerate(node.ancestor_iter()):
            ancestor.append(tmp_node)
            # print(id ,"::::::",tmp_node.index)
            if recom_length < tmp_node.distance_from_tip() :
                attached_node = tmp_node
                attached_id = id
                # print(attached_node.index)
                break

        relocated_nodes = ancestor[attached_id-1]  # relocated node is the adjacent node of recombinant node
        parent.remove_child(node)         # the original recombinant node was removed to reinsert in the other side
        attached_node.remove_child(relocated_nodes) # relocated node was removed to reinsert in to right side
        newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
        newborn.edge_length = attached_node.distance_from_tip() - recom_length
        node.edge_length = recom_length
        newborn.add_child(node)
        relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
        newborn.add_child(relocated_nodes)
        attached_node.add_child(newborn)
        # recombination_tree = tree.as_string(schema="newick")
        recombination_trees.append(tree.as_string(schema="newick"))

    # changeing in topology when recombiantion is larger than tree.max_distance_from_root()
    elif (co_recom + node.edge_length ) >= tree.max_distance_from_root() :
        print(" *********** Stage Three ***********")
        parent.remove_child(node)  # the original recombinant node was removed
        tree.seed_node.add_child(node)
        node.edge_length = co_recom + node.edge_length
        # recombination_tree = tree.as_string(schema="newick")
        recombination_trees.append(tree.as_string(schema="newick"))
    return recombination_trees





def tree_evolver(tree ,node ,nu):
    recombination_trees = []
    co_recom = nu/2
    parent = node.parent_node
    grandparent = parent.parent_node
    print("My node is:" , node.index , node.edge_length)
    print("parent::" , parent.index)
    print("grandparent::", grandparent.index)
    print(parent.distance_from_tip())

    # topology does not change in this case:
    if (co_recom + node.edge_length) < (grandparent.distance_from_tip()):
        print(" *********** Stage one ***********")
        node.edge.length = node.edge.length + co_recom
        sister = node.sister_nodes()
        sister[0].edge.length = sister[0].edge.length + co_recom
        parent.edge.length = parent.edge.length - co_recom
        recombination_trees.append(tree.as_string(schema="newick"))

    # changing in topology to make recombination tree:
    elif ((co_recom + node.edge_length) > grandparent.distance_from_tip())  and ((co_recom + node.edge_length) < tree.max_distance_from_root()):
        print(" *********** Stage Two ***********")
        ancestor = []
        recom_length = co_recom + node.edge_length
        for id,tmp_node in enumerate(node.ancestor_iter()):
            ancestor.append(tmp_node)
            # print(id ,"::::::",tmp_node.index)
            if recom_length < tmp_node.distance_from_tip() :
                attached_node = tmp_node
                attached_id = id
                # print(attached_node.index)
                break

        relocated_nodes = ancestor[attached_id-1]  # relocated node is the adjacent node of recombinant node
        parent.remove_child(node)         # the original recombinant node was removed to reinsert in the other side
        attached_node.remove_child(relocated_nodes) # relocated node was removed to reinsert in to right side
        newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
        newborn.edge_length = attached_node.distance_from_tip() - recom_length
        node.edge_length = recom_length
        newborn.add_child(node)
        relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
        newborn.add_child(relocated_nodes)
        attached_node.add_child(newborn)
        recombination_trees.append(tree.as_string(schema="newick"))

    # changeing in topology when recombiantion is larger than tree.max_distance_from_root()
    elif (co_recom + node.edge_length ) >= tree.max_distance_from_root() :
        print(" *********** Stage Three ***********")
        parent.remove_child(node)  # the original recombinant node was removed
        tree.seed_node.add_child(node)
        node.edge_length = co_recom + node.edge_length
        recombination_trees.append(tree.as_string(schema="newick"))
    return recombination_trees



taxon = tree.taxon_namespace
print(taxon)
pdm = tree.phylogenetic_distance_matrix()
mrca = pdm.mrca(taxon[0], taxon[2])
print(mrca.index)


node0 = tree.find_node_with_taxon_label(label ="0")
# re0 = tree_evolver(tree, node0 , .22 )
# print(re0)
re1 = tree_evolver2(tree, node0 , .6 )
print(re1)
# print(tree.as_ascii_plot(plot_metric='length'))
# print(tree.as_string(schema="newick"))

# tree = Tree.get_from_path(tree_path, 'newick')
# myPhylo.set_index(tree,column[0])
#
# node1 = tree.find_node_with_taxon_label(label ="4")
# re1 = tree_evolver(tree, node1 , .01)
# print(re1)
# print(tree.as_ascii_plot(plot_metric='length'))
#
# tree = Tree.get_from_path(tree_path, 'newick')
# myPhylo.set_index(tree,column[0])
#
#
# re2 = tree_evolver(tree, mrca , 0.8 )
# print(re2)
# print(tree.as_ascii_plot(plot_metric='length'))




































