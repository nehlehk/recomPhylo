import hmmlearn.base
import numpy as np
import logging
import hmmlearn._utils
import myPhylo
import dendropy



def compute_logprob_phylo_old(X, tree, nodes, model):
    recom_trees = myPhylo.make_recombination_trees(tree, 5, nodes)
    n, dim = X.shape
    result = np.zeros((n, len(recom_trees)))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        myPhylo.reroot_tree(state_tree, nodes)
        children = state_tree.seed_node.child_nodes()
        for site_id, partial in enumerate(X):
            p = np.zeros(4)
            p = np.dot(model.p_matrix(children[0].edge_length), partial[0:4])
            for i in range(1, len(children)):
                p *= np.dot(model.p_matrix(children[i].edge_length), partial[i * 4:(i + 1) * 4])
            result[site_id, tree_id] = sum(p)
    return result


def compute_logprob_phylo(X, recom_trees, model):
    n, dim = X.shape
    result = np.zeros((n, len(recom_trees)))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        children = state_tree.seed_node.child_nodes()
        for site_id, partial in enumerate(X):
            p = np.zeros(4)
            p = np.dot(model.p_matrix(children[0].edge_length), partial[0:4])
            for i in range(1, len(children)):
                p *= np.dot(model.p_matrix(children[i].edge_length), partial[i * 4:(i + 1) * 4])
            # result[site_id, tree_id] = sum(p)
            # print(p)
            site_l = np.dot(p, model.get_pi())
            result[site_id, tree_id] = np.log(site_l)
    return result


_log = logging.getLogger(__name__)



class phyloLL_HMM(hmmlearn.base._BaseHMM):
    def __init__(self, n_components, trees, model):
        super().__init__(n_components)
        self.trees = trees
        self.model = model

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
        logprob : array, shape (n_samples, n_components)
            Log probability of each sample in ``X`` for each of the
            model states.
        """

        return compute_logprob_phylo(X, self.trees, self.model)

#     ==========================================================================
