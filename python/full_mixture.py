from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as la
import phyloHMM
import xml.etree.ElementTree as ET




# ==============================================   input  ==============================================================
tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/BaciSim/1/RAxML_bestTree.tree'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/BaciSim/1/wholegenome.fasta"), schema="fasta")


pi = [0.2184,0.2606,0.3265,0.1946]
rates = [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ]
GTR_sample = myPhylo.GTR_model(rates,pi)

column = myPhylo.get_DNA_fromAlignment(alignment)
dna = column[0]
myPhylo.set_index(tree,alignment)
tips = len(dna)

tips_num = len(alignment)
alignment_len = alignment.sequence_size


print("Original tree")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot(show_internal_node_labels = True))

taxon = tree.taxon_namespace
nu = 0.03
# ==============================================   methods  ============================================================
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

        left = np.dot(sqrtPiInv, eigvec)
        right = np.dot(eigvec_inv, sqrtPi)

        p = np.dot(left, np.dot(np.diag(np.exp(eigval * br_length)), right))


        return p


def give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3


def set_index(tree, dna):
    sequence_count = len(dna)

    for node in tree.postorder_node_iter():
        node.index = -1
        node.annotations.add_bound_attribute("index")

    s = sequence_count
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            node.index = s
            node.label = str(node.index)
            s += 1
        else:
            for idx, name in enumerate(dna):
                if str(name) == str(node.taxon):
                    node.index = idx
                    node.label = str(node.index)
                    break


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


def set_tips_partial(tree, alignment):
    partial = np.zeros(((alignment_len, tips_num, 4)))
    for tip in range(tips_num):
      for site in range(alignment_len):
        dna = column[site]
        i = give_index(str(dna[tip]))
        partial[site][tip][i] = 1
    return partial


def computelikelihood_mixture(tree, alignment, tip_partial, model):
    alignment_len = alignment.sequence_size
    tips = len(dna)
    partial = np.zeros(((alignment_len, (2 * tips) - 1, 4)))
    partial[:, 0:tips, :] = tip_partial
    persite_ll = []
    for site in range(alignment_len):
        for node in tree.postorder_node_iter():
            if not node.is_leaf():
                children = node.child_nodes()
                partial[site, node.index] = np.dot(model.p_matrix(children[0].edge_length),
                                                   partial[site, children[0].index])
                for i in range(1, len(children)):
                    partial[site, node.index] *= np.dot(model.p_matrix(children[i].edge_length),
                                                        partial[site, children[i].index])
        p = np.dot(partial[site, tree.seed_node.index], model.get_pi())
        persite_ll.append(np.log(p))

    return persite_ll, partial



def make_hmm_input_mixture(tree, alignment, tip_partial, model):
    sitell, partial = computelikelihood_mixture(tree, alignment, tip_partial, model)
    children = tree.seed_node.child_nodes()
    children_count = len(children)
    x = np.zeros((alignment.sequence_size, children_count * 4))
    for id, child in enumerate(children):
        print(child.index)
        x[:, (id * 4):((id + 1) * 4)] = partial[:, child.index, :]
    return x



def update_mixture_partial(alignment,tree,node,tipdata,posterior,node_order):
  column = get_DNA_fromAlignment(alignment)

  for site in range(alignment_len):
    dna = column[site]
    my_number = give_index(dna[node.index])
    rho = give_rho(node,posterior,site,tips_num,node_order)
    for i in range(4):
      if i == my_number:
        tipdata[site,node.index,i] = 1
      else:
        tipdata[site,node.index,i] = rho

  return tipdata


def tree_evolver_rerooted(tree ,node ,nu):
    co_recom = nu/2
    if (node.edge_length is None):
       node.edge.length = 0
    node.edge.length = node.edge.length + co_recom
    recombination_tree = tree.as_string(schema="newick")

    return recombination_tree

def give_rho(node,recom_prob,site,tips_num,node_order):
  parent = node.parent_node
  if parent == tree.seed_node:
    myindex = parent.index -1
  else:
    myindex = parent.index
  # node_prob = recom_prob[myindex - tips_num][site]
  # node_prob = recom_prob[site]
  # rho = 1 - node_prob[0]
  rho = recom_prob[site][node_order]

  return rho


def give_taxon(tree,index):
  for node in tree.postorder_node_iter():
    if node.index == index:
        return str(node.taxon)




mytree = []
posterior = []
hiddenStates = []
score = []
figure = {}
tipdata = set_tips_partial(tree, alignment)
for id_tree, target_node in enumerate(tree.postorder_internal_node_iter(exclude_seed_node=True)):
    print(target_node.index)
    recombination_trees = []
    child_order = []

    mytree.append(Tree.get_from_path(tree_path, 'newick'))
    set_index(mytree[id_tree], alignment)

    # ----------- Step 1 : Make input for hmm ------------------------------------------------------
    # --------------  Stetp 1.1 : re-root the tree based on the target node where the target node is each internal node of the tree.

    mytree[id_tree].reroot_at_node(target_node, update_bipartitions=False, suppress_unifurcations=True)
    recombination_trees.append(mytree[id_tree].as_string(schema="newick"))

    # --------------  Step 1.2: Calculate X based on this re-rooted tree

    X = make_hmm_input_mixture(mytree[id_tree], alignment, tipdata, GTR_sample)
    # Y = myPhylo.make_hmm_input(mytree[id_tree], alignment, GTR_sample)
    print(X[4500])

    # ----------- Step 2: make 3 recombination trees -----------------------------------------------
    temptree = {}
    recom_child_order = []

    for id, child in enumerate(target_node.child_node_iter()):
        recom_child_order.append(child.index)

    for id, child in enumerate(target_node.child_node_iter()):
        # print(child.index)
        temptree["tree{}".format(id)] = Tree.get_from_path(tree_path, 'newick')
        set_index(temptree["tree{}".format(id)], alignment)

        filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
        target_node_temp = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
        temptree["tree{}".format(id)].reroot_at_node(target_node_temp, update_bipartitions=False,
                                                     suppress_unifurcations=True)

        kids =temptree["tree{}".format(id)].seed_node.child_nodes()
        for kid in kids:
            recom_child_order.append(kid.index)
        filter_fn = lambda n: hasattr(n, 'index') and n.index == child.index
        recombined_node = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
        recombination_trees.append(tree_evolver_rerooted(temptree["tree{}".format(id)], recombined_node, nu))
        child_order.append(recombined_node.index)

    print(child_order)
    print(recom_child_order)
    # ----------- Step 3: Call phyloHMM ----------------------------------------------------------
    model = phyloHMM.phyloLL_HMM(n_components=4, trees=recombination_trees, model=GTR_sample , child_order=child_order , recom_child_order = recom_child_order)
    model.startprob_ = np.array([0.6, 0.13, 0.13, 0.14])
    model.transmat_ = np.array([[0.9997, 0.0001, 0.0001, 0.0001],
                              [0.001, 0.997, 0.001, 0.001],
                              [0.001, 0.001, 0.997, 0.001],
                              [0.001, 0.001, 0.001, 0.997]])


    p = model.predict_proba(X)
    print("posterior:")
    print(p[4500])
    # posterior = model.predict_proba(X)
    # hiddenStates = model.predict(X)
    # score = model.score(X)

    posterior.append(p)
    # hiddenStates.append(model.predict(X))
    # score.append(model.score(X))

    # figure["plot{}".format(id_tree)] = plt.figure(figsize=(15, 8))
    # ax = figure["plot{}".format(id_tree)].add_subplot(2, 1, 1)
    # ax.set_title("Internal Node" + str(target_node.index) + " -- log probability of the most likely state is  " + str(score))
    # ax.plot(hiddenStates)
    # ax.set_ylabel("Clonal - Recombination State")
    #
    # ax2 = figure["plot{}".format(id_tree)].add_subplot(2, 1, 2)
    # ax2.plot(posterior[:, 0], label="ClonalFrame")
    # ax2.plot(posterior[:, 1], label="Recombination A ")
    # ax2.plot(posterior[:, 2], label="Recombination B ")
    # ax2.plot(posterior[:, 3], label="Recombination C ")
    # ax2.set_ylabel("posterior probability for each state")
    # ax2.legend(loc=1, bbox_to_anchor=(1.13, 1.1))

    tree_updatePartial = Tree.get_from_path(tree_path, 'newick')
    set_index(tree_updatePartial, alignment)
    filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
    target_node_partial = tree_updatePartial.find_node(filter_fn=filter_fn)
    for id, child in enumerate(target_node_partial.child_node_iter()):
        if child.is_leaf():
            order = child_order.index(child.index)
            # order = id
            print("my beloved child:", child.index, child.taxon, "order:", order + 1)
            new_partial = update_mixture_partial(alignment, tree_updatePartial, child, tipdata, p ,order+1)

    # print(new_partial[4500])


# tmp = tipdata.transpose(1, 0, 2)
# print(tmp.shape)
#
# my_xml = ET.parse('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/externalRecom/externalRecomTemplate.xml')
# root = my_xml.getroot()
# data = root.find("data")
#
# for i in range(tmp.shape[0]):
#     x = ''
#     c = ET.Element("sequence")
#     c.set("taxon" , give_taxon(tree,i))
#     c.set("uncertain" , "true")
#     for j in range(tmp.shape[1]):
#       x = x + str(repr(tmp[i,j,:]))[7:-2] + ';'
#     c.text = '\n' + x +'\n'
#     data.insert(i,c)
#     c.tail = "\n"
#
# my_xml.write('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/externalRecom/externalRecomPartial.xml' ,encoding="utf-8", xml_declaration=True)


plt.show()