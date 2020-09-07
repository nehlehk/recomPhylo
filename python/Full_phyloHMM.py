from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np
import phyloHMM
import matplotlib.pyplot as plt


# ==============================================   input  ==============================================================
tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/exampledataset/exampledataset_RAxML_bestTree'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/exampledataset/wholegenome.fasta"), schema="fasta")


pi = [0.2184,0.2606,0.3265,0.1946]
rates = [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ]
GTR_sample = myPhylo.GTR_model(rates,pi)

column = myPhylo.get_DNA_fromAlignment(alignment)
dna = column[0]
myPhylo.set_index(tree,dna)


taxon = tree.taxon_namespace
nu = 0.5
# ============================================  methods ================================================================



# ==============================================   input  ==============================================================
# print("Original tree:::::::::::::::")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot())


mytree = []
posterior = []
hiddenStates = []
score = []
for id_tree, target_node in enumerate(tree.postorder_internal_node_iter(exclude_seed_node=True)):
        print(target_node.index)
        # print(id_tree)
        recombination_trees = []
        mytree.append(Tree.get_from_path(tree_path, 'newick'))
        myPhylo.set_index(mytree[id_tree], dna)
        # ----------- Step 1 : Make input for hmm ------------------------------------------------------
        # --------------  Stetp 1.1 : re-root the tree based on the target node where the target node is each internal node of the tree.
        mytree[id_tree].reroot_at_node(target_node, update_bipartitions=False ,suppress_unifurcations = True)
        recombination_trees.append( mytree[id_tree].as_string(schema="newick"))
        # --------------  Step 1.2: Calculate X based on this re-rooted tree
        X = myPhylo.make_hmm_input(mytree[id_tree], alignment, GTR_sample)
        # print(X)

        # ----------- Step 2: make 3 recombination trees -----------------------------------------------
        temptree = {}
        for id, child in enumerate(target_node.child_node_iter()):
                # print(child.index)
                temptree["tree{}".format(id)] = Tree.get_from_path(tree_path, 'newick')
                myPhylo.set_index(temptree["tree{}".format(id)], dna)

                filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
                target_node_temp = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
                temptree["tree{}".format(id)].reroot_at_node(target_node_temp, update_bipartitions=False,suppress_unifurcations=True)

                filter_fn = lambda n: hasattr(n, 'index') and n.index == child.index
                recombined_node = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
                recombination_trees.append(myPhylo.tree_evolver_rerooted(temptree["tree{}".format(id)],recombined_node,nu))


        # ----------- Step 3: Call phyloHMM -----------------------------------------------------
        model = phyloHMM.phyloLL_HMM(n_components=4, trees=recombination_trees,  model=GTR_sample)
        model.startprob_ = np.array([0.85, 0.02, 0.12, 0.01])
        model.transmat_ = np.array( [[0.9999, 0.00002, 0.00006, 0.00002],
                                     [0.0007, 0.999, 0.0002, 0.0001],
                                     [0.0008, 0.0001, 0.999, 0.0001],
                                     [0.0008, 0.0001, 0.0001, 0.999]])

        posterior = model.predict_proba(X)
        hiddenStates = model.predict(X)
        score = model.score(X)

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

# ======================================================================================================================
print("posterior:::::::::::::::::")
print(posterior)

print("hiddenStates::::::::::::::::")
print(hiddenStates)

print("score:::::::::::::::::::::")
print(score)






















# i wrote to this part of code to change the topology!
#     recombination_trees = []
#     temptree = {}
#     target_kids = find_kids_index(target_node)
#     for id, child in enumerate(target_kids):
#         print(child)
#         temptree["tree{}".format(id)] = Tree.get_from_path(tree_path, 'newick')
#         myPhylo.set_index(temptree["tree{}".format(id)], column[0])
#         filter_fn = lambda n: hasattr(n, 'index') and n.index == child
#         recombined_node = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
#         position = find_node_position(recombined_node,target_kids)
#         print(position)
#         tmp = myPhylo.tree_evolver(temptree["tree{}".format(id)],recombined_node,.24 , position)
#         print(tmp)

# for id, child in enumerate(target_node.child_node_iter()):
#     print(child.index)
#     temptree["tree{}".format(id)] = Tree.get_from_path(tree_path, 'newick')
#     myPhylo.set_index(temptree["tree{}".format(id)], column[0])
#     i = child.index
#     filter_fn = lambda n: hasattr(n, 'index') and n.index == i
#     recombined_node = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
#     tmp = myPhylo.tree_evolver(temptree["tree{}".format(id)],recombined_node,.24)
#     print(tmp)