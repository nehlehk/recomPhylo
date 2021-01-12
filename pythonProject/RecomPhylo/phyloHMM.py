import numpy as np
import numpy.linalg as la
import pandas as pd
import matplotlib.pyplot as plt
import hmmlearn.base
import logging
import hmmlearn._utils
from scipy.integrate import quad
from dendropy import Tree, DnaCharacterMatrix
import dendropy
from hmmlearn import hmm
import string
import xml.etree.ElementTree as ET
from dendropy.calculate import treecompare
from sklearn.metrics import mean_squared_error
import math
import argparse




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
# **********************************************************************************************************************
class phyloLL_HMM(hmmlearn.base._BaseHMM):
    def __init__(self, n_components, trees, model ,child_order,recom_child_order):
        super().__init__(n_components)
        self.trees = trees
        self.model = model
        self.child_order = child_order
        self.recom_child_order = recom_child_order

    def _init(self, X, lengths):

        """Initializes model parameters prior to fitting.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Feature matrix of individual samples.
            n_samples == number of alignment sites
            n_features == 12: 4 site partials for each of 3 neighbour nodes

        lengths : array-like of integers, shape (n_sequences, )
            Lengths of the individual sequences in ``X``. The sum of
            these should be ``n_samples``.
        """
        init = 1. / self.n_components
        if 's' in self.init_params or not hasattr(self, "startprob_"):
            self.startprob_ = np.full(self.n_components, init)
        if 't' in self.init_params or not hasattr(self, "transmat_"):
            self.transmat_ = np.full((self.n_components, self.n_components), init)
        n_fit_scalars_per_param = self._get_n_fit_scalars_per_param()
        n_fit_scalars = sum(n_fit_scalars_per_param[p] for p in self.params)
        if X.size < n_fit_scalars:
            _log.warning("Fitting a model with {} free scalar parameters with "
                         "only {} data points will result in a degenerate "
                         "solution.".format(n_fit_scalars, X.size))

    #     ==========================================================================
    def _check(self):
        """Validates model parameters prior to fitting.

        Raises
        ------

        ValueError
            If any of the parameters are invalid, e.g. if :attr:`startprob_`
            don't sum to 1.
        """
        self.startprob_ = np.asarray(self.startprob_)
        if len(self.startprob_) != self.n_components:
            raise ValueError("startprob_ must have length n_components")
        if not np.allclose(self.startprob_.sum(), 1.0):
            raise ValueError("startprob_ must sum to 1.0 (got {:.4f})"
                             .format(self.startprob_.sum()))

        self.transmat_ = np.asarray(self.transmat_)
        if self.transmat_.shape != (self.n_components, self.n_components):
            raise ValueError(
                "transmat_ must have shape (n_components, n_components)")
        if not np.allclose(self.transmat_.sum(axis=1), 1.0):
            raise ValueError("rows of transmat_ must sum to 1.0 (got {})"
                             .format(self.transmat_.sum(axis=1)))

    #     ==========================================================================
    def _compute_log_likelihood(self, X):
        """Computes per-component log probability under the model.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Feature matrix of individual samples.

        Returns
        -------
            Log probability of each sample in ``X`` for each of the
            model states.
        """

        return compute_logprob_phylo(X, self.trees, self.model, self.child_order, self.recom_child_order)

# **********************************************************************************************************************
def give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3
# **********************************************************************************************************************
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
# **********************************************************************************************************************
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
# **********************************************************************************************************************
def set_tips_partial(tree, alignment):
    partial = np.zeros(((alignment_len, tips_num, 4)))
    for tip in range(tips_num):
      for site in range(alignment_len):
        dna = column[site]
        i = give_index(str(dna[tip]))
        partial[site][tip][i] = 1
    return partial
# **********************************************************************************************************************
def computelikelihood_mixture(tree, alignment ,tip_partial, model):
    alignment_len = alignment.sequence_size
    tips = len(dna)
    partial = np.zeros(((alignment_len,(2 * tips) -1, 4)))
    partial[:,0:tips,:] = tip_partial
    persite_ll = np.zeros(alignment_len)


    column = get_DNA_fromAlignment(alignment)

    uniqueCol = list(set(column))
    for u in range(len(uniqueCol)):
      indexes = []
      indexes = [id for id, x in enumerate(column) if x == uniqueCol[u]]
      site = indexes[0]
      for node in tree.postorder_node_iter():
          if not node.is_leaf():
              children = node.child_nodes()
              partial[site,node.index] = np.dot(model.p_matrix(children[0].edge_length), partial[site,children[0].index])
              for i in range(1, len(children)):
                  partial[site, node.index] *= np.dot(model.p_matrix(children[i].edge_length),partial[site,children[i].index])
      partial[indexes,:,:] = partial[site,:,:]
      p = np.dot(partial[site,tree.seed_node.index] , model.get_pi())
      persite_ll[indexes] = np.log(p)


    return persite_ll, partial
# **********************************************************************************************************************
def make_hmm_input_mixture(tree, alignment, tip_partial, model):
    sitell, partial = computelikelihood_mixture(tree, alignment, tip_partial, model)
    children = tree.seed_node.child_nodes()
    children_count = len(children)
    x = np.zeros((alignment.sequence_size, children_count * 4))
    for id, child in enumerate(children):
        # print(child.index)
        x[:, (id * 4):((id + 1) * 4)] = partial[:, child.index, :]
    return x
# **********************************************************************************************************************
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
# **********************************************************************************************************************
def give_rho(node,recom_prob,site,tips_num,node_order):
  parent = node.parent_node
  if parent == tree.seed_node:
    myindex = parent.index -1
  else:
    myindex = parent.index
  rho = recom_prob[site][node_order]
  return rho
# **********************************************************************************************************************
def tree_evolver_rerooted(tree ,node ,nu):
    co_recom = nu/2
    if (node.edge_length is None):
       node.edge.length = 0
    node.edge.length = node.edge.length + co_recom
    recombination_tree = tree.as_string(schema="newick")

    return recombination_tree
# **********************************************************************************************************************
def give_taxon(tree,index):
    for node in tree.postorder_node_iter():
        if int(node.index) == index:
            return int(node.taxon.label)
# **********************************************************************************************************************
def compute_logprob_phylo_fast(X, recom_trees, model):
    n, dim = X.shape
    result = np.zeros((n, len(recom_trees)))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        children = state_tree.seed_node.child_nodes()

        uniqueX = np.unique(X, axis=0)
        for uniq_id, partial in enumerate(uniqueX):
            indexes = []
            indexes = [id for id, x in enumerate(X) if (x == uniqueX[uniq_id]).all()]
            p = np.zeros(4)
            p = np.dot(model.p_matrix(children[0].edge_length), partial[0:4])
            for i in range(1, len(children)):
                p *= np.dot(model.p_matrix(children[i].edge_length), partial[i * 4:(i + 1) * 4])
            site_l = np.dot(p, model.get_pi())
            result[indexes, tree_id] = np.log(site_l)
    return result
# **********************************************************************************************************************
def compute_logprob_phylo(X, recom_trees, model,child_order,recom_child_order):
    n, dim = X.shape
    result = np.zeros((n, len(recom_trees)))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        children = state_tree.seed_node.child_nodes()
        for site_id, partial in enumerate(X):
            order = child_order.index(recom_child_order[tree_id * len(children)])
            p = np.zeros(4)
            p = np.dot(model.p_matrix(children[0].edge_length), partial[order * 4:(order + 1) * 4])
            for i in range(1, len(children)):
                order = child_order.index(recom_child_order[(tree_id* len(children)) + i])
                p *= np.dot(model.p_matrix(children[i].edge_length), partial[order * 4:(order + 1) * 4])
            site_l = np.dot(p, model.get_pi())
            result[site_id, tree_id] = np.log(site_l)
    return result
# **********************************************************************************************************************
def phylohmm(tree,alignment,nu):
    mytree = []
    posterior = []
    hiddenStates = []
    score = []
    tipdata = set_tips_partial(tree, alignment)
    for id_tree, target_node in enumerate(tree.postorder_internal_node_iter(exclude_seed_node=True)):
        # print(target_node.index)
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


        # ----------- Step 2: make recombination trees -----------------------------------------------
        temptree = {}
        recom_child_order = []

        for id, child in enumerate(target_node.child_node_iter()):  # to keep the order of clonal tree children
            recom_child_order.append(child.index)

        for id, child in enumerate(target_node.child_node_iter()):
            # print(child.index)
            temptree["tree{}".format(id)] = Tree.get_from_path(tree_path, 'newick')
            set_index(temptree["tree{}".format(id)], alignment)

            filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
            target_node_temp = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
            temptree["tree{}".format(id)].reroot_at_node(target_node_temp, update_bipartitions=False,suppress_unifurcations=True)

            kids = temptree["tree{}".format(id)].seed_node.child_nodes()  # to keep the order of recombination trees children
            for kid in kids:
                recom_child_order.append(kid.index)

            filter_fn = lambda n: hasattr(n, 'index') and n.index == child.index
            recombined_node = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
            recombination_trees.append(tree_evolver_rerooted(temptree["tree{}".format(id)], recombined_node, nu))
            child_order.append(recombined_node.index)


        # ----------- Step 3: Call phyloHMM -----------------------------------------------------

        model = phyloLL_HMM(n_components=4, trees=recombination_trees, model=GTR_sample, child_order=child_order,
                            recom_child_order=recom_child_order)
        model.startprob_ = np.array([0.6, 0.13, 0.13, 0.14])
        model.transmat_ = np.array([[0.9997, 0.0001, 0.0001, 0.0001],
                                    [0.001, 0.997, 0.001, 0.001],
                                    [0.001, 0.001, 0.997, 0.001],
                                    [0.001, 0.001, 0.001, 0.997]])

        p = model.predict_proba(X)

        posterior.append(p)
        hiddenStates.append(model.predict(X))
        score.append(model.score(X))

        tree_updatePartial = Tree.get_from_path(tree_path, 'newick')
        set_index(tree_updatePartial, alignment)
        filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
        target_node_partial = tree_updatePartial.find_node(filter_fn=filter_fn)
        for id, child in enumerate(target_node_partial.child_node_iter()):
            # print(child.index)
            if child.is_leaf():
                order = child_order.index(child.index)
                # print("my beloved child:", child.index , child.taxon , "order:" , order+1)
                update_mixture_partial(alignment, tree_updatePartial, child, tipdata, p, order + 1)

    return tipdata,posterior,hiddenStates,score
# **********************************************************************************************************************
def ranges(nums):
    nums = sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    return list(zip(edges, edges))
# **********************************************************************************************************************
def recom_ranges(tipdata,threshold):
    my_tipdata = tipdata.transpose(1, 0, 2)
    recom_index = []
    for i in range(my_tipdata.shape[1]):
        for j in range(10):
            if ((my_tipdata[j, i, 0] > threshold) and (my_tipdata[j, i, 0] < 1.0)) or ((my_tipdata[j, i, 1] > threshold) and (my_tipdata[j, i, 1] < 1.0)):
                recom_index.append(i)

    return ranges(recom_index)
# **********************************************************************************************************************
def recom_resultFig(tree,tipdata,threshold):
    my_tipdata = tipdata.transpose(1, 0, 2)
    output = np.zeros((alignment_len, tips_num))
    for i in range(my_tipdata.shape[1]):
        for j in range(my_tipdata.shape[0]):
            if ((my_tipdata[j, i, 0] > threshold) and (my_tipdata[j, i, 0] < 1.0)) or ((my_tipdata[j, i, 1] > threshold) and (my_tipdata[j, i, 1] < 1.0)):
                output[i, j] = 1

    fig = plt.figure(figsize=(15, 3))
    taxa = output.shape[1]
    color = ['red', 'green', 'purple', 'blue', 'black']
    my_taxon = []
    for i in range(taxa):
        my_taxon.append(give_taxon(tree, int(i)))

    t = np.sort(my_taxon)
    for i in range(taxa):
        ax = fig.add_subplot(taxa, 1, i + 1)
        id = my_taxon.index(t[i])
        ax.plot(output[:, id], label=t[i], color=color[i % 5])
        ax.legend(loc=1, bbox_to_anchor=(1.03, 1.4))
        ax.set_frame_on(False)
        ax.axis('off')
    ax.axis('on')
    ax.set_yticklabels([])
    plt.savefig("PhyloHMM_Recombination.jpeg")
# **********************************************************************************************************************
def internal_plot(posterior,hiddenStates,score):
    for i in range(len(posterior)):
        poster = posterior[i]
        hidden = hiddenStates[i]
        sc = score[i]

        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(2, 1, 1)
        ax.set_title("Hidden Markov Models - ClonalFrame and Recombination -- log probability of the most likely state is  " + str(sc))
        ax.plot(hidden)
        ax.set_ylabel("Clonal - Recombination State")
        ax2 = fig.add_subplot(2, 1, 2)
        ax2.plot(poster[:, 0], label="ClonalFrame")
        ax2.plot(poster[:, 1], label="Recombination A ")
        ax2.plot(poster[:, 2], label="Recombination B ")
        ax2.plot(poster[:, 3], label="Recombination C ")
        ax2.set_ylabel("posterior probability for each state")
        ax2.legend(loc=1, bbox_to_anchor=(1.13, 1.1))
        plt.savefig("posterior"+str(i)+".jpeg")
# **********************************************************************************************************************
def make_beast_xml(tipdata,tree,xml_path):
    my_tipdata = tipdata.transpose(1, 0, 2)
    my_xml = ET.parse(xml_path)
    root = my_xml.getroot()
    data = root.find("data")

    for i in range(my_tipdata.shape[0]):
        x = ''
        c = ET.Element("sequence")
        c.set("taxon" , str(give_taxon(tree,i)))
        c.set("uncertain" , "true")
        for j in range(my_tipdata.shape[1]):
          x = x + str(repr(my_tipdata[i,j,:]))[7:-2] + ';'
        c.text = '\n' + x +'\n'
        data.insert(i,c)
        c.tail = "\n"

    my_xml.write('RecomPartial.xml' ,encoding="utf-8", xml_declaration=True)
# **********************************************************************************************************************



if __name__ == "__main__":


    # tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/BaciSim/2/RAxML_bestTree.tree'
    # tree = Tree.get_from_path(tree_path, 'newick')
    # alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/BaciSim/2/wholegenome.fasta"),schema="fasta")
    # xml_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/BaciSim/2/Beast/GTR_template.xml'

    # tree_path = '/home/nehleh/Documents/work/55/b103191f338ba6dcad98e21ad62b23/RAxML_bestTree.tree'
    # tree = Tree.get_from_path(tree_path, 'newick')
    # alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/work/1f/f6627cde194e3110eee2327f16428c/wholegenome.fasta"),schema="fasta")


    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-t', "--treeFile", type=str, required= True, help='tree')
    parser.add_argument('-a', "--alignmentFile", type=str, required= True , help='fasta file')
    parser.add_argument('-nu', "--nuHmm", type=float,default=0.03,help='nuHmm')
    parser.add_argument('-f', "--frequencies", type=list, default= [0.2184,0.2606,0.3265,0.1946],help='frequencies')
    parser.add_argument('-r', "--rates", type=list, default= [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ], help='rates')
    parser.add_argument('-x', "--xmlFile", type=str, required=True, help='xmlFile')
    args = parser.parse_args()

    tree_path = args.treeFile
    genomefile = args.alignmentFile
    xml_path = args.xmlFile
    pi = args.frequencies
    rates = args.rates
    nu = args.nuHmm




    tree = Tree.get_from_path(tree_path, 'newick')
    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    GTR_sample = GTR_model(rates, pi)
    column = get_DNA_fromAlignment(alignment)
    dna = column[0]
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size

    set_index(tree, alignment)
    tipdata,posterior,hiddenStates,score = phylohmm(tree, alignment, nu)

    internal_plot(posterior, hiddenStates, score)

    tree = Tree.get_from_path(tree_path, 'newick')
    set_index(tree, alignment)
    recom_resultFig(tree,tipdata, 0.3)

    make_beast_xml(tipdata,tree,xml_path)

