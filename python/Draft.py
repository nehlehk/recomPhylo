from dendropy import Tree, DnaCharacterMatrix
import dendropy
import myPhylo
import numpy as np
import phyloHMM


tree = Tree.get_from_path('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/tree_6taxa.tree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/phyloHMM/sample_6taxa.fasta"), schema="fasta")


pi = [0.2184,0.2606,0.3265,0.1946]
rates = [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ]


GTR_sample = myPhylo.GTR_model(rates,pi)


# sitell , partial =myPhylo.wholeAlignmentLikelihood(tree,alignment, GTR_sample)


print("Before Re-rooting:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())



# mrca = reroot_tree(tree, [0,1])
# print("After Re-rooting:")
# print(tree.as_string(schema='newick'))
# print(tree.as_ascii_plot())



x =  myPhylo.make_hmm_input(tree,alignment,GTR_sample,[0,2])
# print(x[0])
# print(x[7])
print(x)
# print(x.shape)




# recom_trees = make_recombination_trees(tree,0.05,[0,2])
# print(recom_trees)































