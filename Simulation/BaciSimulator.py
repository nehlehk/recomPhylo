import numpy as np
import dendropy
import math
from dendropy.simulate import treesim
import random
from dendropy import Tree
import matplotlib.pyplot as plt
from operator import itemgetter, attrgetter
import math


tips_number = 10
alignment_len = 5000
recom_len = 500
threshold_len = 200
recom_rate = .05
max_tMRCA= 0.01
nu_ex = 0.25
nu_in = 0.1


taxon_list= []
for i in range(tips_number):
  taxon_list.append(str(i))
taxa = dendropy.TaxonNamespace(taxon_list)
tree = treesim.pure_kingman_tree(taxon_namespace=taxa,pop_size=3)
# print(tree.as_string(schema="newick"))
# print(tree.as_ascii_plot())
nodes_number = len(tree.nodes())


t = tree.max_distance_from_root()
normal_co = t/ max_tMRCA
# print(normal_co)
tree.scale_edges(1/normal_co)
clonal_tree = tree.as_string(schema="newick")
print(tree.as_string(schema="newick"))
print(tree.as_ascii_plot())
clonal_tree = clonal_tree.replace('\n',"")


# Poisson( tree.sum() * rel_recomb_rate_per_site * alignment_length)
recom_num =  np.random.poisson(tree.length() * recom_rate * alignment_len)

print(recom_num)

# ----------------------------------------------------------------------------------------------------------------------
def ex_recom_maker(tree ,node ,nu):
    rand_nu = np.random.normal(nu,0.01)
    co_recom = rand_nu/2

    if co_recom > tree.max_distance_from_root():
      if (node.edge_length is None):
        node.edge.length = 0
      if node.is_leaf():
        # changeing in topology when recombiantion is larger than tree.max_distance_from_root()
        if (co_recom + node.edge_length ) >= tree.max_distance_from_root() :
            new_tree = dendropy.Tree(taxon_namespace=taxa)
            external_len = co_recom + node.edge_length - tree.max_distance_from_root()
            # tree2 = tree.extract_tree_without_taxa_labels(node.taxon.label)
            tree2 = tree.extract_tree_without_taxa_labels(node.label)
            other_nodes = tree2.seed_node
            new_tree.seed_node.add_child(other_nodes)
            new_tree.seed_node.add_child(node)
            other_nodes.edge_length = external_len
            node.edge_length = co_recom + node.edge_length
            return new_tree.as_string(schema="newick")
      else:
         # changeing in topology when recombiantion is larger than tree.max_distance_from_root()
        if (co_recom + node.edge_length ) >= tree.max_distance_from_root() :
            new_tree = dendropy.Tree(taxon_namespace=taxa)
            # external_len = co_recom + node.edge_length - tree.max_distance_from_root()
            # external_len = co_recom + node.edge_length + tree.max_distance_from_root()
            tree.prune_subtree(node)
            other_nodes = tree.seed_node
            new_tree.seed_node.add_child(other_nodes)
            new_tree.seed_node.add_child(node)
            # other_nodes.edge_length = external_len
            node.edge_length = co_recom + node.edge_length
            other_nodes.edge_length = node.distance_from_tip() + node.edge_length - other_nodes.distance_from_tip()

            return new_tree.as_string(schema="newick")

    else:
      print("This nu value can not make an external recombination! Please set the nu value bigger than tMRCA.")
      return tree
# ----------------------------------------------------------------------------------------------------------------------
def set_index(tree):
    for node in tree.postorder_node_iter():
        node.index = -1
        node.annotations.add_bound_attribute("index")

    s = len(tree.leaf_nodes())
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            node.index = s
            node.label = str(node.index)
            s += 1
        else:
            node.index = node.taxon.label
            node.label = str(node.index)
# ----------------------------------------------------------------------------------------------------------------------
def give_descendents(tree, node_index):
    if node_index >= tips_number:
        internal_recom_node = tree.find_node_with_label(str(node_index))
        children = internal_recom_node.child_nodes()
        # print(children)
        for n in range(len(children)):
            r_node = int(children[n].index)
            if r_node >= tips_number:
                give_descendents(tree, r_node)
            else:
                desc.append(r_node)

    return desc
# ----------------------------------------------------------------------------------------------------------------------
def find_recom_tree(index,site):
  for i in range(recom_num):
    if (all_data[i,0] == str(float(index))) and ((site >= float(all_data[i,1])) and (site <= float(all_data[i,3]))):
      return all_data[i,4]
# ----------------------------------------------------------------------------------------------------------------------
def give_equivalent_node(recomtree):
  e = []
  for edge in recomtree.postorder_edge_iter():
    if edge.length is None:
      edge.length = 0
    e.append(edge.length)

  m = max(e)
  for node in recomtree.postorder_node_iter():
    if node.edge_length == m:
      return node.label,m
# ----------------------------------------------------------------------------------------------------------------------
def merge_trees_new(recomtrees, recomnodes):
    equ = np.zeros((len(recomtrees), 4))
    for treeid in range(len(recomtrees)):
        rtree = Tree.get_from_string(recomtrees[treeid], schema='newick')
        set_index(rtree)
        equ[treeid, 0] = recomnodes[treeid]
        equ[treeid, 1:3] = give_equivalent_node(rtree)
        equ[treeid, 3] = treeid

    # print(equ)
    s_equ = equ[equ[:, 2].argsort()[::-1]]

    # s_equ = sorted(equ, key=itemgetter(1), reverse=False)
    print(s_equ)

    clonaltree = Tree.get_from_string(clonal_tree, schema='newick')
    set_index(clonaltree)

    print(clonaltree.as_ascii_plot())
    print(clonaltree.as_string(schema="newick"))
    for node in clonaltree.postorder_node_iter():
        print(node.index)

    for i in range(len(recomtrees)):
        # print(clonaltree.as_ascii_plot())
        # print(clonaltree.as_string(schema="newick"))

        temptree = Tree.get_from_string(recomtrees[int(s_equ[i][3])], schema='newick')
        set_index(temptree)
        # print(temptree.as_ascii_plot())
        # print(temptree.as_string(schema="newick"))
        prunenode = clonaltree.find_node_with_label(str(int(s_equ[i][0])))
        print("prunenode:")
        print(prunenode)
        print(prunenode.edge_length)

        node = temptree.find_node_with_label(str(int(s_equ[i][1])))
        print("node:")
        print(node)
        print(node.edge_length)

        parent = prunenode.parent_node

        # changing in topology to make recombination tree:
        # if ((node.edge_length) > parent.edge_length + prunenode.edge_length) and ((node.edge_length) < clonaltree.max_distance_from_root()):
        if ((node.edge_length) < clonaltree.max_distance_from_root()):
            # print(" *********** Stage Two ***********")
            ancestor = []
            for id, tmp_node in enumerate(prunenode.ancestor_iter()):
                ancestor.append(tmp_node)
                # print(id ,"::::::",tmp_node.index)
                if node.edge_length < tmp_node.distance_from_tip():
                    attached_node = tmp_node
                    attached_id = id
                    # print(attached_node.index)
                    break

            relocated_nodes = ancestor[attached_id - 1]  # relocated node is the adjacent node of recombinant node
            # parent.remove_child(prunenode)         # the original recombinant node was removed to reinsert in the other side
            clonaltree.prune_subtree(prunenode)

            # print(clonaltree.as_ascii_plot())
            # print(clonaltree.as_string(schema="newick"))

            attached_node.remove_child(relocated_nodes)  # relocated node was removed to reinsert in to right side
            newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
            # newborn.edge_length = attached_node.distance_from_tip() - node.edge_length
            # print("-----------------------------------:",node.distance_from_tip())
            # newborn.add_child(node)
            newborn.add_child(prunenode)
            prunenode.edge_length = node.edge_length
            newborn.edge_length = clonaltree.max_distance_from_root() - (node.edge_length + node.distance_from_tip())
            relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
            newborn.add_child(relocated_nodes)
            attached_node.add_child(newborn)
            # print(clonaltree.as_ascii_plot())
            # print(clonaltree.as_string(schema="newick"))

        # changeing in topology when recombiantion is larger than tree.max_distance_from_root()
        if node.edge_length > clonaltree.max_distance_from_root():
            # print(" *********** Stage Three ***********")
            new_tree = dendropy.Tree(taxon_namespace=taxa)
            if node.is_leaf():
                print("is leaf")
                tree2 = clonaltree.extract_tree_without_taxa_labels(prunenode.label)
                other_nodes = tree2.seed_node
                new_tree.seed_node.add_child(other_nodes)
                new_tree.seed_node.add_child(node)
                other_nodes.edge_length = node.edge_length - clonaltree.max_distance_from_root()
                clonaltree = new_tree
                print(clonaltree.as_ascii_plot())
                print(clonaltree.as_string(schema="newick"))
            else:
                print("is internal")
                # clonaltree.prune_subtree(node)
                clonaltree.prune_subtree(prunenode)
                other_nodes = clonaltree.seed_node
                new_tree.seed_node.add_child(other_nodes)
                new_tree.seed_node.add_child(node)
                other_nodes.edge_length = node.edge_length + node.distance_from_tip() - other_nodes.distance_from_tip()
                clonaltree = new_tree
                for node in clonaltree.postorder_node_iter():
                    print(node.index)
                print(clonaltree.as_ascii_plot())
                print(clonaltree.as_string(schema="newick"))

    return clonaltree.as_string(schema="newick")
    # print(clonaltree.as_ascii_plot())
    # print(clonaltree.as_string(schema="newick"))
# ----------------------------------------------------------------------------------------------------------------------


set_index(tree)
node_labels = []
node_weight = []
for node in tree.postorder_node_iter():
  if not node==tree.seed_node:
    node_labels.append(node.index)
    node_weight.append(node.edge_length /tree.length())
ideal_gap = int(alignment_len/recom_num)
info = np.zeros((recom_num,4))
my_trees = []
recom_info = []
for recom_id in range(int(recom_num)):
  starting_falg = False
  tree = Tree.get_from_string(clonal_tree,schema='newick')
  set_index(tree)
  while not starting_falg:
    random_tip = np.random.choice(node_labels,1,node_weight)
    start_pos = random.randint(ideal_gap*recom_id, ideal_gap * (recom_id+1) )
    r_len = random.randint(threshold_len, recom_len + threshold_len)
    if (start_pos + r_len <= alignment_len) and (int(random_tip)!= info[recom_id-1,0]):
      info[recom_id,0] = int(random_tip)
      info[recom_id,1] = start_pos
      info[recom_id,2] = r_len
      info[recom_id,3] = start_pos + r_len
      recom_node = tree.find_node_with_label(random_tip)
      recom_tree= ex_recom_maker(tree,recom_node,nu_ex) # make external recombination
      my_trees.append(recom_tree)
      starting_falg = True

my_trees = np.array(my_trees).reshape(-1,1)
all_data = np.append(info,my_trees, 1)




output = np.zeros((alignment_len, nodes_number))
for i in range(nodes_number):
  for j in range(recom_num):
    if all_data[j,0] == str(float(i)):
      s = int(float(all_data[j,1]))
      e = int(float(all_data[j,3]))
      output[s:e,i] = 1






myoutput = output.transpose()
site = 0
sum = 1
over_recom_pos = 0
while site < alignment_len -1:
  clonal = 0
  recom = 0
  overlap = 0
  recom_node = []
  overcount = 0
  sum_over = 0
  result_tree = ''
  over_trees = []
  while (np.count_nonzero(myoutput[:,site]) == 0) and (site < alignment_len-1) :
    clonal = clonal + 1
    site = site + 1
  if clonal > 0:
    print("clonal:" ,clonal , clonal_tree)
  while (np.count_nonzero(myoutput[:,site]) == 1)  and (site < alignment_len-1):
      recom_node = np.where(myoutput[:,site]==1)[0]
      index = np.where(all_data[:,1] == str(float(site)))[0]
      recom = recom + 1
      site = site + 1
      # if len(index)== 0:
      #   index = np.where(all_data[:,1] == str(float(over_recom_pos)))[0]
      if len(index)== 0:
        index = np.where(all_data[:,3] == str(float(site)))[0]
      if len(index) > 0:
        result_tree = all_data[int(index),4]
  if recom > 0:
    print("recom:",recom ,"recome_node:",recom_node , result_tree)
  while (np.count_nonzero(myoutput[:,site]) > 1) and (site < alignment_len-1):
    overlap = 0
    over_trees = []
    over_nodes = []
    mergetree = ''
    over_count = np.count_nonzero(myoutput[:,site])
    while  (np.count_nonzero(myoutput[:,site]) == over_count):
       over_nodes = np.where(myoutput[:,site]==1)[0]
       overlap = overlap + 1
       site = site + 1
    if len(over_nodes) > 0:
      for t in range(len(over_nodes)):
        over_trees.append(find_recom_tree(over_nodes[t],site))

      print(over_nodes,"--------",over_trees)
      # mergetree = merge_trees_new(over_trees,over_nodes)
      print("overlap:",overlap,"over_nodes:",over_nodes , mergetree)
  # over_recom_pos = site - overlap


  sum = sum + clonal + recom + overlap
print(sum)


fig = plt.figure(figsize=(15,7))
color=['red', 'green' ,'purple', 'blue','black']

for i in range(nodes_number):
  ax = fig.add_subplot(nodes_number,1,i+1)
  ax.plot(output[:,i] ,label= i ,color = color[i%5])
  ax.legend(loc = 1 , bbox_to_anchor=(1.03, 1.4))
  ax.set_frame_on(False)
  ax.axis('off')

ax.axis('on')
ax.set_yticklabels([])
plt.show()