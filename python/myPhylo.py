import numpy
import numpy.linalg as la


class GTR_model:
    def __init__(self, rates, pi):
        self.rates = rates
        self.pi = pi

    #     ========================================================================
    def p_matrix(self , br_length):
        p = numpy.zeros((4, 4))

        mu = 0
        freq = numpy.zeros((4, 4))
        q = numpy.zeros((4, 4))
        sqrtPi = numpy.zeros((4, 4))
        sqrtPiInv = numpy.zeros((4, 4))
        exchang = numpy.zeros((4, 4))
        s = numpy.zeros((4, 4))
        fun = numpy.zeros(1)
        a, b, c, d, e = self.rates
        f = 1

        freq = numpy.diag(self.pi)
        sqrtPi = numpy.diag(numpy.sqrt(self.pi))
        sqrtPiInv = numpy.diag(1.0 / numpy.sqrt(self.pi))
        mu = 1 / (2 * ((a * self.pi[0] * self.pi[1]) + (b * self.pi[0] * self.pi[2]) + (c * self.pi[0] * self.pi[3]) + (d * self.pi[1] * self.pi[2]) + (
                e * self.pi[1] * self.pi[3]) + (self.pi[2] * self.pi[3])))
        exchang[0][1] = exchang[1][0] = a
        exchang[0][2] = exchang[2][0] = b
        exchang[0][3] = exchang[3][0] = c
        exchang[1][2] = exchang[2][1] = d
        exchang[1][3] = exchang[3][1] = e
        exchang[2][3] = exchang[3][2] = f

        q = numpy.multiply(numpy.dot(exchang, freq), mu)

        for i in range(4):
            q[i][i] = -sum(q[i][0:4])

        s = numpy.dot(sqrtPi, numpy.dot(q, sqrtPiInv))

        eigval, eigvec = la.eig(s)
        eigvec_inv = la.inv(eigvec)

        left = numpy.dot(sqrtPi, eigvec)
        right = numpy.dot(eigvec_inv, sqrtPiInv)

        p = numpy.dot(left, numpy.dot(numpy.diag(numpy.exp(eigval * br_length)), right))

        return p
    #=============================================================================
def computelikelihood(tree, dna , model):

    tips = len(dna)
    partial = numpy.zeros(((2 * tips) -1, 4))
    set_index(tree,dna)

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
            partial[node.index] = numpy.dot(model.p_matrix(children[0].edge_length), partial[children[0].index])
            for i in range(1, len(children)):
                # print("child index:" ,children[i].index , "edge_length  : ", children[i].edge_length)
                partial[node.index] *= numpy.dot(model.p_matrix(children[i].edge_length),
                                                 partial[children[i].index])



    # print(partial)
    # ll_partial = numpy.zeros(2* tips)
    # for i in range(tree.seed_node.index, tips -1, -1):
    #     ll_partial[i] = round(numpy.log(numpy.mean(partial[i]) ), 7)
    #
    #
    # persite_ll = ll_partial[tree.seed_node.index]
    persite_ll = round(numpy.log(numpy.mean(partial[tree.seed_node.index]) ), 7)

    return persite_ll, partial

#     ========================================================================
def partial_likelihoods_to_target_node(tree,dna,model):

    tips = len(dna)
    partial = numpy.zeros((2 * tips, 4))
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
            partial[node.index] = numpy.dot(model.p_matrix(children[0].edge_length), partial[children[0].index])
            for i in range(1, len(children)):
                partial[node.index] *= numpy.dot(model.p_matrix(children[i].edge_length),
                                                 partial[children[i].index])



    for i in range(partial.shape[0]):
        print(numpy.mean(partial[i,]))

    ll_partial = numpy.zeros(2)
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            children = node.child_nodes()
            for i in range(0, len(children)):
                # print(children[i].index)
                ll_partial[i] = round(numpy.log(numpy.mean(partial[children[i].index])), 7)

    # print(ll_partial)
    return ll_partial

#     ========================================================================
def expectedLL(tree, dna, co_clonal, co_recom, model):

    tips = len(dna)

    partial = numpy.zeros((2 * tips, 4))
    exp_clonal = numpy.zeros((2 * tips, 4))
    exp_recom = numpy.zeros((2 * tips, 4))

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
            partial[node.index] = numpy.dot(model.p_matrix(children[0].edge_length), partial[children[0].index])
            exp_clonal[node.index] = numpy.dot(model.p_matrix(children[0].edge_length * co_clonal),
                                               partial[children[0].index])
            exp_recom[node.index] = numpy.dot(model.p_matrix(children[0].edge_length * co_recom),
                                              partial[children[0].index])
            for i in range(1, len(children)):
                partial[node.index] *= numpy.dot(model.p_matrix(children[i].edge_length),
                                                 partial[children[i].index])
                exp_clonal[node.index] *= numpy.dot(model.p_matrix(children[i].edge_length * co_clonal),
                                                    exp_clonal[children[i].index])
                exp_recom[node.index] *= numpy.dot(model.p_matrix(children[i].edge_length * co_recom),
                                                   exp_recom[children[i].index])

    expected_clonal_ll = numpy.zeros(tips - 1)
    expected_recombination_ll = numpy.zeros(tips - 1)
    p_index = 0
    for par in range(tips + 1, tree.seed_node.index + 1, 1):
        expected_clonal_ll[p_index] = round(numpy.log(numpy.mean(exp_clonal[par])), 7)
        expected_recombination_ll[p_index] = round(numpy.log(numpy.mean(exp_recom[par])), 7)
        p_index += 1


    return expected_clonal_ll, expected_recombination_ll
#     =======================================================================================
def wholeAlignmentExpLL(tree, alignment, internl_node=-1, co_clonal=1, co_recom=1):

    tips = len(alignment)
    alignment_len = alignment.sequence_size

    clonal_vector = numpy.zeros((alignment_len, tips -1))
    recom_vector = numpy.zeros((alignment_len, tips -1))

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
    LL_root = numpy.zeros(alignment_len)
    LL_partial = numpy.zeros(((alignment_len , (2*tips)-1 , 4)))

    column = get_DNA_fromAlignment(alignment)

    uniqueCol = list(set(column))
    for i in range(len(uniqueCol)):
        indexes = [id for id, x in enumerate(column) if x == uniqueCol[i]]
        result = computelikelihood(tree, uniqueCol[i] , model)
        LL_root[indexes] = result[0]
        LL_partial[indexes,:] = result[1]

    return LL_root , LL_partial

#     =======================================================================================
def make_recombination_trees1(tree, co_recom):
    recombination_trees = []
    recombination_trees.append(tree.as_string(schema="newick"))
    for node in tree.postorder_node_iter():
        if node.edge.length is None:
            node.edge.length = 0
        # print(node.edge.length)
        if node.edge.length > 0:
            node.edge.length = node.edge.length * co_recom
            recombination_trees.append(tree.as_string(schema="newick"))
            node.edge.length = node.edge.length * 1 / co_recom

    return recombination_trees

#     =======================================================================================
def make_recombination_trees(tree, dna,targetnode,co_recom):
    set_index(tree,dna)
    recombination_trees = []
    recombination_trees.append(tree.as_string(schema="newick"))
    for node in tree.postorder_node_iter():
        if node.edge.length is None:
            node.edge.length = 0
        # print(node.edge.length)
        if (node.edge.length > 0) and not(node.parent_node.index == targetnode)  :
            node.edge.length = node.edge.length * co_recom
            recombination_trees.append(tree.as_string(schema="newick"))
            node.edge.length = node.edge.length * (1 / co_recom)

    return recombination_trees



