from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np
import phyloHMM
import matplotlib.pyplot as plt



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


taxon = tree.taxon_namespace
nu = 0.03

print(alignment.sequence_size)


print("Original tree:::::::::::::::")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot(show_internal_node_labels = True))








mytree = []
recombination_trees = []
tipdata = myPhylo.set_tips_partial(tree, alignment)
# print(tipdata[5000])
filter_fn = lambda n: hasattr(n, 'index') and n.index == 16
target_node = tree.find_node(filter_fn=filter_fn)
# print(target_node.index)
# ----------- Step 1 : Make input for hmm ------------------------------------------------------
# --------------  Stetp 1.1 : re-root the tree based on the target node where the target node is each internal node of the tree.

tree.reroot_at_node(target_node, update_bipartitions=False ,suppress_unifurcations = True)
recombination_trees.append(tree.as_string(schema="newick"))

# --------------  Step 1.2: Calculate X based on this re-rooted tree

X = myPhylo.make_hmm_input_mixture(tree, alignment, tipdata, GTR_sample)

# Y = myPhylo.make_hmm_input(tree, alignment, GTR_sample)
# X1 = myPhylo.make_hmm_input_new(tree, alignment, GTR_sample)
# mat = myPhylo.make_hmm_input_mathiue(tree, alignment, GTR_sample)




# ----------- Step 2: make 3 recombination trees -----------------------------------------------
temptree = {}
for id, child in enumerate(target_node.child_node_iter()):
        # print(child.index , child.taxon)
        temptree["tree{}".format(id)] = Tree.get_from_path(tree_path, 'newick')
        myPhylo.set_index(temptree["tree{}".format(id)], alignment)
        filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
        target_node_temp = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
        temptree["tree{}".format(id)].reroot_at_node(target_node_temp, update_bipartitions=False,suppress_unifurcations=True)
        filter_fn = lambda n: hasattr(n, 'index') and n.index == child.index
        recombined_node = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
        recombination_trees.append(myPhylo.tree_evolver_rerooted(temptree["tree{}".format(id)],recombined_node,nu))

# print(recombination_trees)
# ----------- Step 3: Call phyloHMM -----------------------------------------------------

model = phyloHMM.phyloLL_HMM(n_components=4, trees=recombination_trees,  model=GTR_sample)
# model.startprob_ = np.array([0.85, 0.05, 0.05, 0.05])
model.startprob_ = np.array([0.94, 0.02, 0.02, 0.02])
model.transmat_ = np.array([[0.9999999, 0.00000003, 0.00000003, 0.00000004],
                            [0.00003, 0.9999, 0.00003, 0.00004],
                            [0.00004, 0.00003, 0.9999, 0.00003],
                            [0.00003, 0.00004, 0.00003, 0.9999]])

posterior = model.predict_proba(X)
hiddenStates = model.predict(X)
score = model.score(X)



# print(posterior[1])

tree_updatePartial = Tree.get_from_path(tree_path, 'newick')
myPhylo.set_index(tree_updatePartial, alignment)
filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
target_node_partial = tree_updatePartial.find_node(filter_fn=filter_fn)
for id, child in enumerate(target_node_partial.child_node_iter()):
    if child.is_leaf():
        # print("my beloved child:", child.index , child.taxon)
        new_partial = myPhylo.update_mixture_partial(alignment, tree_updatePartial, child, tipdata, posterior)

# print(tipdata[5000])


fig = plt.figure(figsize=(15, 8))

ax = fig.add_subplot(2, 1, 1)
ax.set_title("Hidden Markov Models - ClonalFrame and Recombination -- log probability of the most likely state is  " + str(score))
ax.plot(hiddenStates)
ax.set_ylabel("Clonal - Recombination State")

ax2 = fig.add_subplot(2, 1, 2)
ax2.plot(posterior[:, 0], label="ClonalFrame")
ax2.plot(posterior[:, 1], label="Recombination A ")
ax2.plot(posterior[:, 2], label="Recombination B ")
ax2.plot(posterior[:, 3], label="Recombination C ")
ax2.set_ylabel("posterior probability for each state")
ax2.legend(loc=1, bbox_to_anchor=(1.13, 1.1))



plt.show()