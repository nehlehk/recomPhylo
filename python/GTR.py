from dendropy import Tree, DnaCharacterMatrix
import numpy
import math
import dendropy
import numpy.linalg as la


class GTR:


    def give_index(c):
        if c == "A":
            return 0
        elif c == "C":
            return 1
        elif c == "G":
            return 2
        elif c == "T":
            return 3




    def p_matrix(br_length):
        p = numpy.zeros((4, 4))

        mu = 0
        freq = numpy.zeros((4, 4))
        q = numpy.zeros((4, 4))
        sqrtPi = numpy.zeros((4, 4))
        sqrtPiInv = numpy.zeros((4, 4))
        exchang = numpy.zeros((4, 4))
        s = numpy.zeros((4, 4))
        fun = numpy.zeros(1)
        a, b, c, d, e = rates
        f = 1

        freq = numpy.diag(pi)
        sqrtPi = numpy.diag(numpy.sqrt(pi))
        sqrtPiInv = numpy.diag(1.0 / numpy.sqrt(pi))
        mu = 1 / (2 * ((a * pi[0] * pi[1]) + (b * pi[0] * pi[2]) + (c * pi[0] * pi[3]) + (d * pi[1] * pi[2]) + (
                e * pi[1] * pi[3]) + (pi[2] * pi[3])))
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




    def computelikelihood(tree, dna):

        partial = numpy.zeros((2 * tips, 4))

        for node in tree.postorder_node_iter():
            node.index = -1
            node.annotations.add_bound_attribute("index")

        s = tips + 1
        for id, node in enumerate(tree.postorder_node_iter()):
            if not node.is_leaf():
                node.index = s
                s += 1
            else:
                for idx, name in enumerate(dna):
                    if idx + 1 == int(node.taxon.label):
                        node.index = idx + 1
                        break
        pos = 0
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                i = give_index(str(dna[pos]))
                pos += 1
                # i = give_index(dna[node.index-1])
                partial[node.index][i] = 1
            else:
                children = node.child_nodes()
                partial[node.index] = numpy.dot(p_matrix(children[0].edge_length, model), partial[children[0].index])
                for i in range(1, len(children)):
                    partial[node.index] *= numpy.dot(p_matrix(children[i].edge_length, model),
                                                     partial[children[i].index])


        ll_partial = numpy.zeros(tips - 1)
        p_index = 0
        for par in range(tips + 1, tree.seed_node.index + 2, 1):
            temp = 0
            temp = numpy.sum(partial[par]) * 0.25
            # ll_partial.append(round(numpy.log(temp),7))
            ll_partial[p_index] = round(numpy.log(temp), 7)
            p_index += 1

        # ll = numpy.sum(partial[12]) * 0.25

        # ll = numpy.sum(partial[tree.seed_node.index]) * 0.25

        # print("likelihood = {} and log-likelihood = {} ".format(round(ll,7) , round(numpy.log(ll),7)))
        # return round(numpy.log(ll),7)
        # return 1
        return ll_partial