import numpy as np
import numpy.linalg as la
import dendropy
from dendropy import Tree


class GTR_model:
    def __init__(self, rates, pi):
        self.rates = rates
        self.pi = pi
    #     ========================================================================
    def get_pi(self):
        return self.pi
    #     ========================================================================
    def p_matrix(self , br_length):
        p = np.zeros((4, 4))

        mu = 0
        freq = np.zeros((4, 4))
        q = np.zeros((4, 4))
        sqrtPi = np.zeros((4, 4))
        sqrtPiInv = np.zeros((4, 4))
        exchang = np.zeros((4, 4))
        s = np.zeros((4, 4))
        fun = np.zeros(1)
        a, b, c, d, e = self.rates
        f = 1

        freq = np.diag(self.pi)
        sqrtPi = np.diag(np.sqrt(self.pi))
        sqrtPiInv = np.diag(1.0 / np.sqrt(self.pi))
        mu = 1 / (2 * ((a * self.pi[0] * self.pi[1]) + (b * self.pi[0] * self.pi[2]) + (c * self.pi[0] * self.pi[3]) + (d * self.pi[1] * self.pi[2]) + (
                e * self.pi[1] * self.pi[3]) + (self.pi[2] * self.pi[3])))
        exchang[0][1] = exchang[1][0] = a
        exchang[0][2] = exchang[2][0] = b
        exchang[0][3] = exchang[3][0] = c
        exchang[1][2] = exchang[2][1] = d
        exchang[1][3] = exchang[3][1] = e
        exchang[2][3] = exchang[3][2] = f


        q = np.multiply(np.dot(exchang, freq), mu)

        for i in range(4):
            q[i][i] = -sum(q[i][0:4])


        s = np.dot(sqrtPi, np.dot(q, sqrtPiInv))


        eigval, eigvec = la.eig(s)
        eigvec_inv = la.inv(eigvec)


        # left = np.dot(sqrtPi, eigvec)
        # right = np.dot(eigvec_inv, sqrtPiInv)

        left = np.dot(sqrtPiInv, eigvec)
        right = np.dot(eigvec_inv, sqrtPi)

        p = np.dot(left, np.dot(np.diag(np.exp(eigval * br_length)), right))


        return p
    #=============================================================================
def computelikelihood(tree, dna , model):

    tips = len(dna)
    partial = np.zeros(((2 * tips) -1, 4))
    # set_index(tree,dna)

    pos = 0
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            i = give_index(str(dna[pos]))
            pos += 1
            partial[node.index][i] = 1
        else:
            # print("node.index: ",node.index)
            children = node.child_nodes()
            # print("child index:" ,children[0].index , "  edge_length0 : ",children[0].edge_length)
            partial[node.index] = np.dot(model.p_matrix(children[0].edge_length), partial[children[0].index])
            for i in range(1, len(children)):
                # print("child index:" ,children[i].index , "edge_length  : ", children[i].edge_length)
                partial[node.index] *= np.dot(model.p_matrix(children[i].edge_length),partial[children[i].index])




    p = np.dot(partial[tree.seed_node.index] , model.get_pi())
    persite_ll = np.log(p)

    return persite_ll, partial
#     ========================================================================
def partial_likelihoods_to_target_node(tree,dna,model):

    tips = len(dna)
    partial = np.zeros((2 * tips, 4))
    set_index(tree, dna)

    pos = 0
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            i = give_index(str(dna[pos]))
            pos += 1
            partial[node.index][i] = 1 * node.edge.length
        else:
            # print(node.index)
            children = node.child_nodes()
            partial[node.index] = np.dot(model.p_matrix(children[0].edge_length), partial[children[0].index])
            for i in range(1, len(children)):
                partial[node.index] *= np.dot(model.p_matrix(children[i].edge_length),
                                                 partial[children[i].index])



    for i in range(partial.shape[0]):
        print(np.mean(partial[i,]))

    ll_partial = np.zeros(2)
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            children = node.child_nodes()
            for i in range(0, len(children)):
                # print(children[i].index)
                ll_partial[i] = round(np.log(np.mean(partial[children[i].index])), 7)

    # print(ll_partial)
    return ll_partial

#     ========================================================================
def expectedLL(tree, dna, co_clonal, co_recom, model):

    tips = len(dna)

    partial = np.zeros((2 * tips, 4))
    exp_clonal = np.zeros((2 * tips, 4))
    exp_recom = np.zeros((2 * tips, 4))

    set_index(tree,dna)

    pos = 0
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            i = give_index(str(dna[pos]))
            pos += 1
            partial[node.index][i] = 1
            exp_clonal[node.index][i] = 1
            exp_recom[node.index][i] = 1
        else:
            children = node.child_nodes()
            partial[node.index] = np.dot(model.p_matrix(children[0].edge_length), partial[children[0].index])
            exp_clonal[node.index] = np.dot(model.p_matrix(children[0].edge_length * co_clonal),
                                               partial[children[0].index])
            exp_recom[node.index] = np.dot(model.p_matrix(children[0].edge_length * co_recom),
                                              partial[children[0].index])
            for i in range(1, len(children)):
                partial[node.index] *= np.dot(model.p_matrix(children[i].edge_length),
                                                 partial[children[i].index])
                exp_clonal[node.index] *= np.dot(model.p_matrix(children[i].edge_length * co_clonal),
                                                    exp_clonal[children[i].index])
                exp_recom[node.index] *= np.dot(model.p_matrix(children[i].edge_length * co_recom),
                                                   exp_recom[children[i].index])

    expected_clonal_ll = np.zeros(tips - 1)
    expected_recombination_ll = np.zeros(tips - 1)
    p_index = 0
    for par in range(tips + 1, tree.seed_node.index + 1, 1):
        expected_clonal_ll[p_index] = round(np.log(np.mean(exp_clonal[par])), 7)
        expected_recombination_ll[p_index] = round(np.log(np.mean(exp_recom[par])), 7)
        p_index += 1


    return expected_clonal_ll, expected_recombination_ll
#     =======================================================================================
def wholeAlignmentExpLL(tree, alignment, internl_node=-1, co_clonal=1, co_recom=1):

    tips = len(alignment)
    alignment_len = alignment.sequence_size

    clonal_vector = np.zeros((alignment_len, tips -1))
    recom_vector = np.zeros((alignment_len, tips -1))

    column = get_DNA_fromAlignment(alignment)

    uniqueCol = list(set(column))
    for i in range(len(uniqueCol)):
        indexes = [id for id, x in enumerate(column) if x == uniqueCol[i]]
        result = expectedLL(tree, uniqueCol[i], co_clonal, co_recom)
        clonal_vector[indexes, : ] = result[0]
        recom_vector[indexes, : ] = result[1]

    return clonal_vector, recom_vector
#     =======================================================================================
def give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3
#     ========================================================================
def set_index(tree,dna):
    # result = get_DNA_fromAlignment(alignment)
    # dna = result[0]
    tips = len(dna)
    for node in tree.postorder_node_iter():
        node.index = 0
        node.annotations.add_bound_attribute("index")

    s = tips
    for id, node in enumerate(tree.postorder_node_iter()):
        if not node.is_leaf():
            node.index = s
            s += 1
        else:
            for idx, name in enumerate(dna):
                if idx + 1 == int(node.taxon.label):
                    node.index = idx + 1
                    break

    # for node in tree.postorder_node_iter():
    #     print(node.index)
#     =======================================================================================
def get_DNA_fromAlignment(alignment):
    alignment_len = alignment.sequence_size
    tips = len(alignment)
    column = []
    for l in range(alignment_len):
        col = ""
        for t in range(tips):
            col += str(alignment[t][l])
        column.append(col)

    return column
#     =======================================================================================
def wholeAlignmentLikelihood(tree, alignment , model):
    '''
    :param tree:
    :param alignment:
    :return: persite likelihood and partial likelihood for each site
    '''
    tips = len(alignment)
    alignment_len = alignment.sequence_size
    LL_root = np.zeros(alignment_len)
    LL_partial = np.zeros(((alignment_len , (2*tips)-1 , 4)))

    column = get_DNA_fromAlignment(alignment)

    uniqueCol = list(set(column))
    for i in range(len(uniqueCol)):
        indexes = [id for id, x in enumerate(column) if x == uniqueCol[i]]
        result = computelikelihood(tree, uniqueCol[i] , model)
        LL_root[indexes] = result[0]
        LL_partial[indexes,:] = result[1]

    return LL_root , LL_partial

#     =======================================================================================
def reroot_tree(tree, nodes):
    pdm = tree.phylogenetic_distance_matrix()
    taxon = tree.taxon_namespace
    mrca = pdm.mrca(taxon[nodes[0]], taxon[nodes[1]])
    tree.reroot_at_node(mrca, update_bipartitions=False)
    return mrca
#     =======================================================================================
def make_hmm_input(tree, alignment, model):
    sitell, partial = wholeAlignmentLikelihood(tree, alignment, model)
    children = tree.seed_node.child_nodes()
    children_count = len(children)
    x = np.zeros((alignment.sequence_size, children_count * 4))
    for id, child in enumerate(children):
        # print(child.index)
        x[:, (id * 4):((id + 1) * 4)] = partial[:, child.index, :]
    return x
#     =======================================================================================
def make_recombination_trees_old(tree,co_recom ,params , target_type):
    recombination_trees = []
    recombination_trees.append(tree.as_string(schema="newick"))

    if (target_type == 1) :
        tmp_tree1 = tree.extract_tree()
        mrca = reroot_tree(tmp_tree1, params)
        children = mrca.child_nodes()
        children[0].edge.length = children[0].edge.length + co_recom
        children[1].edge.length = children[1].edge.length + co_recom
        children[2].edge.length = children[2].edge.length - co_recom
        recombination_trees.append(tmp_tree1.as_string(schema="newick"))
        recombination_trees.append(tmp_tree1.as_string(schema="newick"))

        mrca = reroot_tree(tree, params)
        children = mrca.child_nodes()
        children[0].edge.length = children[0].edge.length - co_recom
        children[1].edge.length = children[1].edge.length - co_recom
        children[2].edge.length = children[2].edge.length + co_recom
        recombination_trees.append(tree.as_string(schema="newick"))
    elif (target_type == 2):
        tmp_tree1 = tree.extract_tree()
        mrca = reroot_tree(tmp_tree1, params)
        children = mrca.child_nodes()
        tmp_tree1.prune_taxa(children[0].taxon)


    return recombination_trees
#     =======================================================================================
def tree_evolver(tree ,node ,nu , position):
    recombination_trees = []
    co_recom = nu/2
    parent = node.parent_node
    # grandparent = parent.parent_node
    # print("My node is:" , node.index , node.edge_length)
    # print("parent::" , parent.index)
    # print("grandparent::", grandparent.index)
    # print(parent.distance_from_tip())

    if (node.edge_length is None):
       node.edge.length = 0
    if (parent.edge.length is None):
        parent.edge.length = 0

    # topology does not change in this case:
    if ((co_recom + node.edge_length) < parent.distance_from_tip()) and (position == "descendant"):
        print(" *********** Stage one --- descendant  ***********")
        node.edge.length = node.edge.length + co_recom
        sister = node.sister_nodes()
        sister[0].edge.length = sister[0].edge.length + co_recom
        parent.edge.length = parent.edge.length - co_recom
        recombination_trees.append(tree.as_string(schema="newick"))
    elif ((co_recom + node.edge_length) < parent.distance_from_tip()) and (position == "ancestor"):
        print(" *********** Stage one --- ancestor  ***********")
        node.edge.length = node.edge.length + co_recom

    # changing in topology to make recombination tree:
    elif ((co_recom + node.edge_length) > parent.distance_from_tip())  and ((co_recom + node.edge_length) < tree.max_distance_from_root()):
        # print(" *********** Stage Two ***********")
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
        # print(" *********** Stage Three ***********")
        parent.remove_child(node)  # the original recombinant node was removed
        tree.seed_node.add_child(node)
        node.edge_length = co_recom + node.edge_length
        recombination_trees.append(tree.as_string(schema="newick"))
    return recombination_trees
#     =======================================================================================
def tree_evolver_rerooted(tree ,node ,nu):
    co_recom = nu/2
    if (node.edge_length is None):
       node.edge.length = 0
    node.edge.length = node.edge.length + co_recom
    recombination_tree = tree.as_string(schema="newick")

    return recombination_tree
#     =======================================================================================
def find_kids_index(node):
    kids = []
    for id, child in enumerate(node.child_node_iter()):
        kids.append(child.index)
    return kids
#     =======================================================================================
def find_sibling_index(node):
    s = []
    sibling = node.sibling_nodes()
    for i in range(len(sibling)):
        s.append(sibling[i].index)
    return s
#     =======================================================================================
def find_node_position(node,target_kids):
    sister = find_sibling_index(node)
    if sister[0] in target_kids:
        return "descendant"
    else:
        return "ancestor"
#     =======================================================================================
def make_recombination_trees(tree_path,tree,dna,target_node , nu):
    temptree = {}
    recombination_trees = []
    tree.reroot_at_node(target_node, update_bipartitions=False, suppress_unifurcations=True)
    recombination_trees.append(tree.as_string(schema="newick"))
    for id, child in enumerate(target_node.child_node_iter()):
            temptree["tree{}".format(id)] = Tree.get_from_path(tree_path, 'newick')
            set_index(temptree["tree{}".format(id)], dna)
            temptree["tree{}".format(id)].reroot_at_node(target_node, update_bipartitions=False, suppress_unifurcations=True)
            filter_fn = lambda n: hasattr(n, 'index') and n.index == child.index
            recombined_node = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
            recombination_trees.append(tree_evolver_rerooted(temptree["tree{}".format(id)],recombined_node,nu))
    return recombination_trees
# ============================================================================================
def set_tips_partial(tree, alignment):
    alignment_len = alignment.sequence_size
    tips_num = len(alignment)
    column = get_DNA_fromAlignment(alignment)
    partial = np.zeros(((alignment_len, tips_num, 4)))
    for site in range(alignment_len):
      pos = 0
      for node in tree.postorder_node_iter():
        dna = column[site]
        if node.is_leaf():
          # print(node.index)
          i = give_index(str(dna[pos]))
          pos += 1
          partial[site,node.index,i] = 1
    return partial
# ============================================================================================
def computelikelihood_mixture(tree, alignment, tip_partial, model):
    alignment_len = alignment.sequence_size
    column = get_DNA_fromAlignment(alignment)
    dna = column[0]
    tips = len(dna)
    partial = np.zeros(((alignment_len, (2 * tips) - 1, 4)))
    partial[:, 0:tips, :] = tip_partial
    persite_ll = []
    for site in range(alignment_len):
        for node in tree.postorder_node_iter():
            if not node.is_leaf():
                children = node.child_nodes()
                partial[site, node.index] = np.dot(model.p_matrix(children[0].edge_length),partial[site, children[0].index])
                for i in range(1, len(children)):
                    partial[site, node.index] *= np.dot(model.p_matrix(children[i].edge_length),partial[site, children[i].index])
        p = np.dot(partial[site, tree.seed_node.index], model.get_pi())
        persite_ll.append(np.log(p))

    return persite_ll, partial
# ============================================================================================
def make_hmm_input_mixture(tree, alignment, tip_partial, model):
    sitell, partial = computelikelihood_mixture(tree, alignment, tip_partial, model)
    children = tree.seed_node.child_nodes()
    children_count = len(children)
    x = np.zeros((alignment.sequence_size, children_count * 4))
    for id, child in enumerate(children):
        # print(child.index)
        x[:, (id * 4):((id + 1) * 4)] = partial[:, child.index, :]
    return x
