import dendropy
import numpy as np
import math
import numpy.linalg as la

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
rates = np.ones(5)
pi = [0.25]*4


tree = dendropy.Tree.get_from_string('((1,2)5,3,4)6;', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/0_Research/PhD/Data/test.fasta"), schema="fasta")
partial = np.zeros((4,2*tips))




n = alignment.sequence_size
mu = 0
freq = np.zeros((4, 4))
q = np.zeros((4, 4))
sqrtPi = np.zeros((4, 4))
sqrtPiInv = np.zeros((4, 4))
exchang = np.zeros((4, 4))
s = np.zeros((4, 4))
fun = np.zeros(n)
a, b, c, d, e = rates
f = 1

freq = np.diag(pi)
sqrtPi = np.diag(np.sqrt(pi))
sqrtPiInv = np.diag(1.0 / np.sqrt(pi))
mu = 1 / (2 * ((a * pi[0] * pi[1]) + (b * pi[0] * pi[2]) + (c * pi[0] * pi[3]) + (d * pi[1] * pi[2]) + (
        e * pi[1] * pi[3]) + (pi[2] * pi[3])))
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

left = np.dot(sqrtPi, eigvec)
right = np.dot(eigvec_inv, sqrtPiInv)

p = np.dot(left, np.dot(np.diag(np.exp(eigval * br_length)), right))


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



temp = []
for i in range(4):
    temp.append(partial[i,2*tips -3] * pi[i])
ll = sum(temp)

print("likelihood = {} and log-likelihood = {} ".format(round(ll,7) , round(np.log(ll),7)))