from dendropy import Tree, DnaCharacterMatrix
import numpy as np
import math


def give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3
# =======================================================================================================
dna = 'CGCC'

tips = 4 #tips number
br_length = 0.1
p = np.zeros((4, 4))

# tree = Tree.get_from_path('/home/nehleh/0_Research/PhD/Data/tree.tre', 'newick') (3,4,(1,2)5)6;
tree = Tree.get_from_string('((1,2)5,3,4)6;', 'newick')

# internal_idx = len(tree.taxon_namespace)

# print(internal_idx)

partial = np.zeros((4,2*tips))

for i in range(4):
    for j in range(4):
        p[i][j] = 0.25 - 0.25* math.exp(-4*br_length/3)
    p[i][i] = 0.25 + 0.75*math.exp(-4*br_length/3)
# pprint.pprint(p)


for node in tree.postorder_node_iter():
    node.index = -1
    node.annotations.add_bound_attribute("index")

s = tips + 1
for id,node in enumerate(tree.postorder_node_iter()):
    if not node.is_leaf():
        node.index = s
        s += 1
    else:
        for idx, name in enumerate(dna):
            if idx + 1 == int(node.taxon.label):
                node.index = idx+1
                break


for node in tree.postorder_node_iter():
    if node.is_leaf():
        i = give_index(dna[node.index-1])
        partial[i][node.index-1] = 1
    else:
            for j in range(4):
                sump = []
                for x in node.child_node_iter():
                    z = 0
                    for k in range(4):
                       z  += p[j][k] * partial[k][x.index-1]
                    sump.append(z)
                partial[j][node.index-1] = np.prod(sump)


# print(partial)

ll = sum(partial[:,2*tips -3] * 0.25)

print("likelihood = {} and log-likelihood = {} ".format(round(ll,7) , round(np.log(ll),7)))








