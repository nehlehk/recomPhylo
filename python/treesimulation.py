import dendropy
from dendropy.simulate import treesim
import os

tips_number = 10

temp = []
for i in range(tips_number):
    temp.append(str(i+1))


# taxa = dendropy.TaxonNamespace(temp)
# tree = treesim.birth_death_tree(birth_rate=1.0 , death_rate= 0.5 ,ntax= tips_number , taxon_namespace = taxa)
# tree.print_plot()
# plant = tree.as_string(schema='newick')
# myplant = plant.replace("[&R]", "")
# treefile= "/home/nehleh/0_Research/PhD/Data/simulationdata/tree.tree"
# f = open(treefile, "w")
# f.write(myplant)
# f.close()

treefile = "/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/test-seq-gen/tree.tree"

model = 'GTR'
seqnum = '22000'
frequencies = '0.2184,0.2606,0.3265,0.1946'
rates = '2.0431,0.0821,0,0.067,0,1'
# frequencies = '0.25,0.25,0.25,0.25'
# rates = '1,1,1,1,1,1'
outfile = '/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/test-seq-gen/testseg-gen'
seqgen = '/home/nehleh/0_Research/Software/Seq-Gen-1.3.4/source/seq-gen'
cmd = seqgen + '  -m'+model+ '  -l'+ seqnum + '  -f'+frequencies + '  -r'+rates +'  -of'+ '  <'+treefile+'>  ' + outfile
print(cmd)
os.system(cmd)




