import numpy as np
import dendropy
import math
from dendropy.simulate import treesim
import random
from dendropy import Tree
import matplotlib.pyplot as plt
from operator import itemgetter, attrgetter
import math
import pandas as pd
from collections import Counter


tips_number = 10
alignment_len = 5000
recom_len = 500
threshold_len = 200
recom_rate = .05
max_tMRCA= 0.01
nu_ex = 0.2
# nu_in = 0.1
internal = 1
leaf = 1

taxon_list= []
for i in range(tips_number):
  taxon_list.append(str(i))
taxa = dendropy.TaxonNamespace(taxon_list)


clonal_tree = '(((4:0.001239110283681363,1:0.001239110283681363):0.0018547139845565177,(0:0.0012213226202040935,5:0.0012213226202040935):0.0018725016480337874):0.006906175731762117,(((3:0.0006986298786672515,9:0.0006986298786672515):0.0029400514371517444,6:0.0036386813158189955):0.0024892092060997675,(8:0.0019325752681498172,(2:1.86987151994185e-05,7:1.86987151994185e-05):0.0019138765529503985):0.004195315253768946):0.003872109478081235):0.0;'
recomtrees = ['((((4:0.001239110283681363,1:0.001239110283681363)10:0.0018547139845565177,0:0.003093824268237881)12:0.006906175731762117,(((3:0.0006986298786672515,9:0.0006986298786672515)13:0.0029400514371517444,6:0.0036386813158189955)14:0.0024892092060997675,(8:0.0019325752681498172,(2:1.86987151994185e-05,7:1.86987151994185e-05)15:0.0019138765529503985)16:0.004195315253768946)17:0.003872109478081235)18:0.08750634189181826,5:0.09750634189181825);\n',
 '(((4:0.001239110283681363,1:0.001239110283681363)10:0.008760889716318634,(((3:0.0006986298786672515,9:0.0006986298786672515)13:0.0029400514371517444,6:0.0036386813158189955)14:0.0024892092060997675,(8:0.0019325752681498172,(2:1.86987151994185e-05,7:1.86987151994185e-05)15:0.0019138765529503985)16:0.004195315253768946)17:0.003872109478081235)18:0.09474906003601989,(0:0.0012213226202040935,5:0.0012213226202040935)11:0.1035277374158158);\n']
recomnodes = [5,11]


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
def resolve_parent_nodes(s_equ, tree, temptree):
    desc = []
    n = []
    length = []
    tip_len = []
    for i in range(len(s_equ)):
        my_node = int(s_equ[i][0])
        # print(my_node)
        if my_node >= tips_number:
            res = set()
            desc.append(give_descendents(tree, my_node, res))
            n.append(my_node)
            length.append(s_equ[i][2])
            tip_len.append(temptree.find_node_with_label(str(int(s_equ[i][1]))).distance_from_tip())

    for i in range(len(desc)):
        for j in range(len(s_equ)):
            if (int(s_equ[j][0]) < tips_number) and (int(s_equ[j][0]) in desc[i]):
                tmp = list(desc[i])
                # print(int(s_equ[j][0]) , list(desc[i]) , n[i])
                for k in range(len(tmp)):
                    print(tmp[k], tmp[k], length[i] + tip_len[i], n[i])

    return desc
# ----------------------------------------------------------------------------------------------------------------------
def resolve_conflict(clonaltree, recomtrees, recomnodes):
    conflict = common_nodes(clonaltree, recomnodes)
    print(recomnodes)
    print(conflict)
    if len(conflict) > 0:
        for i in range(len(conflict)):
            # print(conflict[i][0], conflict[i][1])
            node1 = clonaltree.find_node_with_label(str(conflict[i][0]))
            node2 = clonaltree.find_node_with_label(str(conflict[i][1]))
            if node1.parent_node == node2:
                id_main = recomnodes.index(conflict[i][1])
                id_temp = recomnodes.index(conflict[i][0])
                maintree = Tree.get_from_string(recomtrees[id_main], schema='newick')
                set_index(maintree)
                temptree = Tree.get_from_string(recomtrees[id_temp], schema='newick')
                set_index(temptree)

                if int(node1.index) < tips_number:
                    # print("node is taxa")
                    prunenode = maintree.find_node_with_label(str(node1.index))
                    recomnode = temptree.find_node_with_label(str(node1.index))
                else:
                    # print("node is internal")
                    res = set()
                    desc = give_descendents(clonaltree, node1.index, res)
                    prune_index = my_mrca(maintree, list(desc))
                    prunenode = maintree.find_node_with_label(str(prune_index))
                    recom_index = my_mrca(temptree, list(desc))
                    recomnode = temptree.find_node_with_label(str(recom_index))

                # print(maintree.as_ascii_plot())
                # print(maintree.as_string(schema="newick"))
                # print("prunenode:")
                # print(prunenode)
                # print(prunenode.edge_length)
                # print("recomnode:")
                # print(recomnode)
                # print(recomnode.edge_length)
                parent = prunenode.parent_node
                # print("parent")
                # print(parent)

                if ((recomnode.edge_length) > parent.distance_from_tip()) and ((recomnode.edge_length) < maintree.max_distance_from_root()):
                    maintree.prune_subtree(prunenode)  # the original recombinant node was removed to reinsert in the other side

                    res = set()
                    desc = give_descendents(clonaltree, node2.index, res)
                    d = list(desc)
                    random_node = np.random.choice(taxon_list,replace=False)
                    while int(random_node) in d:
                        random_node = np.random.choice(taxon_list,replace=False)
                    candidatenode = maintree.find_node_with_label(str(random_node))

                    attached_node, attached_id, ancestor = find_attached_node(maintree, recomnode, candidatenode)

                    relocated_nodes = ancestor[attached_id - 1]  # relocated node is the adjacent node of recombinant node
                    attached_node.remove_child(relocated_nodes)  # relocated node was removed to reinsert in to right side
                    newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
                    newborn.add_child(recomnode)
                    if attached_node.edge_length is None:
                        attached_node.edge_length = 0
                    newborn.edge_length = maintree.max_distance_from_root() - (recomnode.edge_length + attached_node.edge_length)
                    relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
                    newborn.add_child(relocated_nodes)
                    attached_node.add_child(newborn)



                    # print(maintree.as_string(schema="newick"))
                    # print(maintree.as_ascii_plot())



    return maintree.as_string(schema="newick")
# ----------------------------------------------------------------------------------------------------------------------
def common_nodes(clonaltree,overlap_nodes):
  kids = []
  nodes = []
  s_overlap_nodes = np.sort(overlap_nodes)
  for i in range(len(s_overlap_nodes)):
    nodes.append(s_overlap_nodes[i])
    if s_overlap_nodes[i] >= tips_number:
      desc = set()
      kids.append(give_descendents(clonaltree,s_overlap_nodes[i],desc))
    else:
      kids.append({(s_overlap_nodes[i])})

  subsets = []
  for i in range(len(kids)):
    for j in range(i+1,len(kids)):
      if (kids[i].issubset(kids[j])):
        # print(i,j)
        subsets.append([nodes[i], nodes[j]])

  return subsets
# ----------------------------------------------------------------------------------------------------------------------
def give_descendents(tree,node_index,result):
  if node_index >= tips_number:
    internal_recom_node = tree.find_node_with_label(str(node_index))
    children = internal_recom_node.child_nodes()
    for n in range(len(children)):
      r_node= int(children[n].index)
      if r_node >= tips_number:
        give_descendents(tree,r_node,result)
      else:
        result.add(r_node)
  return result

# ----------------------------------------------------------------------------------------------------------------------
def merge_Recomtrees(recomtrees , recomnodes):

  clonaltree = Tree.get_from_string(clonal_tree,schema='newick')
  set_index(clonaltree)

  # if one one node is the subset of the other nodes
  conflict = common_nodes(clonaltree,recomnodes)

  if (len(conflict) > 0)  and (conflict[0][0] == conflict[0][1]):
        print("Same nodes")
        conflict = []

  if len(conflict) == 0 :
      print("NOOOOOO conflict!")

      desc = []
      equ = np.zeros((len(recomtrees),4))
      for treeid in range(len(recomtrees)):
        rtree= Tree.get_from_string(recomtrees[treeid],schema='newick')
        set_index(rtree)
        equ[treeid,0] = recomnodes[treeid]
        equ[treeid,1:3] = give_equivalent_node(rtree)
        equ[treeid,3] = treeid

      s_equ = equ[equ[:,2].argsort()[::-1]]
      s_equ = give_older_recom(s_equ)
      # print(s_equ)

      maintree = Tree.get_from_string(recomtrees[int(s_equ[0][3])],schema='newick')
      set_index(maintree)

      for i in range(1,len(s_equ)):
        temptree = Tree.get_from_string(recomtrees[int(s_equ[i][3])],schema='newick')
        set_index(temptree)
        node_maintree = int(s_equ[i][0])
        if node_maintree  >= tips_number:
          node_maintree = new_mrca(clonaltree,maintree,node_maintree,int(s_equ[i-1][0]))
        prunenode = maintree.find_node_with_label(str(node_maintree))
        node = temptree.find_node_with_label(str(int(s_equ[i][1])))
        parent = prunenode.parent_node

        attached_node,attached_id,ancestor =  find_attached_node(maintree,node,prunenode)
        regraft_recomnodes(maintree,attached_node,attached_id,prunenode,node,ancestor)

      return maintree.as_string(schema="newick")

  elif len(conflict) > 0 :
    print("CONFLICT!! ")
    for i in range(len(conflict)):
      if conflict[i][0] < tips_number:
        print("one node is taxa!")
      else:
        print("two nodes are internal!")
# ----------------------------------------------------------------------------------------------------------------------
def new_mrca(tree1,tree2,primenode,deletenode):
  a = []
  desc = set()
  d = give_descendents(tree1,primenode,desc)
  rlist = set()
  removelist = give_descendents(tree1,deletenode,rlist)
  if len(removelist)> 0 :
    a = [x for x in d if x not in removelist]
  else:
    a = list(d)
  final_node = my_mrca(tree2,a)
  return final_node

# ----------------------------------------------------------------------------------------------------------------------
def my_mrca(tree,tips):
  pdm = tree.phylogenetic_distance_matrix()
  taxon = tree.taxon_namespace
  node0 = [i for i,x in enumerate(taxon) if x.label==str(tips[0])]
  node1 = [i for i,x in enumerate(taxon) if x.label==str(tips[len(tips)-1])]
  myMrca = pdm.mrca(taxon[node0[0]], taxon[node1[0]])

  return myMrca.index

# ----------------------------------------------------------------------------------------------------------------------

def regraft_recomnodes(tree,attached_node,attached_id,prunenode,recomnode,ancestor):
  relocated_nodes = ancestor[attached_id-1]  # relocated node is the adjacent node of recombinant node
  tree.prune_subtree(prunenode) # the original recombinant node was removed to reinsert in the other side
  attached_node.remove_child(relocated_nodes) # relocated node was removed to reinsert in to right side
  newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
  newborn.add_child(prunenode)
  prunenode.edge_length = recomnode.edge_length
  if attached_node.edge_length is None:
    attached_node.edge_length = 0
  newborn.edge_length = tree.max_distance_from_root() - (recomnode.edge_length + attached_node.edge_length)
  relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
  newborn.add_child(relocated_nodes)
  attached_node.add_child(newborn)
#   --------------------------------------------------------------------------------------------------------------------
def give_equivalent_node(recomtree):
  e = []
  for edge in recomtree.postorder_edge_iter():
    if edge.length is None:
      edge.length = 0
    e.append(edge.length)

  m = np.max(e)
  for node in recomtree.postorder_node_iter():
    if node.edge_length == m and node.is_leaf():
      return node.label,m
    elif node.edge_length == m and node.is_internal():
      return node.label,node.distance_from_tip() + m

# ----------------------------------------------------------------------------------------------------------------------
def find_attached_node(tree,node,prunenode):
  if ((node.edge_length) < tree.max_distance_from_root()) :
    ancestor = []
    for id,tmp_node in enumerate(prunenode.ancestor_iter()):
        ancestor.append(tmp_node)
        if node.edge_length < tmp_node.distance_from_tip() :
            attached_node = tmp_node
            attached_id = id
            break
  return attached_node,attached_id,ancestor

# ----------------------------------------------------------------------------------------------------------------------
def give_older_recom(s_equ):
  count = Counter(a[0] for a in s_equ)
  max = np.max([a[2] for a in s_equ])
  [a for a in s_equ if (count[a[0]] == 1) or (a[2]==max)]
  return [a for a in s_equ if (count[a[0]] == 1) or (a[2]==max)]

# ----------------------------------------------------------------------------------------------------------------------
def recom_maker(tree ,node ,recom_length):
    new_tree = dendropy.Tree(taxon_namespace=taxa)
    parent = node.parent_node
    # changeing in topology when recombiantion is larger than tree.max_distance_from_root()
    if recom_length > tree.max_distance_from_root():
      print("------------------------ recom_length > tree.max_distance_from_root() ------------------------------")
      if (node.edge_length is None):
        node.edge.length = 0
      # new_tree = dendropy.Tree(taxon_namespace=taxa)
      if node.is_leaf():
        external_len = recom_length + node.edge_length - tree.max_distance_from_root()
        tree.prune_subtree(node)
        other_nodes = tree.seed_node
        new_tree.seed_node.add_child(other_nodes)
        new_tree.seed_node.add_child(node)
        other_nodes.edge_length = external_len
        node.edge_length = recom_length + node.edge_length
      else:
        tree.prune_subtree(node)
        other_nodes = tree.seed_node
        new_tree.seed_node.add_child(other_nodes)
        new_tree.seed_node.add_child(node)
        node.edge_length = recom_length + node.edge_length
        other_nodes.edge_length = node.distance_from_tip() + node.edge_length - other_nodes.distance_from_tip()
    #*************************************************************************
    # topology does not change in this case:
    elif ((recom_length + node.edge_length) < parent.distance_from_tip()) and (node.is_leaf):
        print(" *********** ((recom_length + node.edge_length) < parent.distance_from_tip()) and (node.is_leaf)  ***********")
        node.edge.length = node.edge.length + recom_length
        sister = node.sister_nodes()
        sister[0].edge.length = sister[0].edge.length + recom_length
        parent.edge.length = parent.edge.length - recom_length
        new_tree = tree
    elif ((recom_length + node.edge_length) < parent.distance_from_tip()) and (node.is_internal):
        print(" *********** ((recom_length + node.edge_length) < parent.distance_from_tip()) and (node.is_internal)  ***********")
        node.edge.length = node.edge.length + recom_length
        new_tree = tree
    #*************************************************************************
    # changing in topology to make recombination tree:
    elif ((recom_length + node.edge_length) > parent.distance_from_tip())  and ((recom_length + node.edge_length) < tree.max_distance_from_root()):
        print("((recom_length + node.edge_length) > parent.distance_from_tip())  and ((recom_length + node.edge_length) < tree.max_distance_from_root())")
        ancestor = []
        recom_length = recom_length + node.edge_length
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
        new_tree = tree

    #*************************************************************************

    return new_tree


# ----------------------------------------------------------------------------------------------------------------------



initialTree = Tree.get_from_string(clonal_tree,schema='newick')
set_index(initialTree)

merge_Recomtrees(recomtrees,recomnodes)

unconfilcted_str = resolve_conflict(initialTree,recomtrees,recomnodes)
print(unconfilcted_str)

# unconfilcted_tree = Tree.get_from_string(unconfilcted_str, schema='newick')
# set_index(unconfilcted_tree)
#
# tree = Tree.get_from_string(clonal_tree, schema='newick')
# set_index(tree)
#
# res = set()
# desc = give_descendents(tree, 16, res)
# print(desc)
# new_node = my_mrca(unconfilcted_tree, list(desc))
# print(new_node)
#
# res2 = set()
# desc2 = give_descendents(unconfilcted_tree, new_node, res2)
# print(desc2)
#
# mynode = unconfilcted_tree.find_node_with_label(str(new_node))
# result = recom_maker(unconfilcted_tree,mynode,0.10416229)
#
# print(result.as_string(schema="newick"))
# print(result.as_ascii_plot())