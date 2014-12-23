import warnings

from Bio.Phylo import BaseTree, Newick, NewickIO
dir(NewickIO)

#import sys
#import unittest
#from io import StringIO

from Bio import Phylo
#from Bio.Phylo import PhyloXML, NewickIO

def lookup_by_names(tree):
    names = {}
    for clade in tree.get_terminals():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

 
EX_NEWICK = 'spec.nwk'
treeA = Phylo.read(EX_NEWICK, 'newick')
print(treeA)
names = lookup_by_names(treeA)
treeA.prune(names['Bcl'])  #prune() takes an object from the same tree
treeA.count_terminals()
treeA.prune(names['Bha'])  #prune() takes an object from the same tree
treeA.count_terminals()
treeA.prune(names['Bsu'])  #prune() takes an object from the same tree
treeA.count_terminals()
treeA.prune(names['Bpu'])  #prune() takes an object from the same tree
treeA.count_terminals()


help(Phylo.BaseTree) #this is very useful, but lengthy and confusing on screen. 

EX_NEWICK = 'qincodeml/tree.nwk'
tree1 = Phylo.read(EX_NEWICK, 'newick')
tree2 = Phylo.read(EX_NEWICK, 'newick')

#tree2 = tree #This links tree1 and tree2 to the same object
names1 = lookup_by_names(tree1)

print(tree1)
leaf_nodes = tree1.get_terminals()
leaf_nodes[0:3]
len(leaf_nodes)

leaf_nodes2 = tree2.get_terminals()
type(leaf_nodes2[0])   #<class 'Bio.Phylo.Newick.Clade'>

tree2.prune(leaf_nodes2[0])  #prune() takes an object from the same tree
len(tree2.get_terminals())
len(tree1.get_terminals())
tree2.count_terminals()
tree1.count_terminals()

all_clades = tree2.find_clades()
bsu_clades = tree2.find_clades(target="Bsu")
type(bsu_clades) #this is an iterator
# :returns: an iterable through all matching objects, searching
#          depth-first (preorder) by default.

all_clades[bsu_clades]

bsu_node.next #not working? 

tree2.prune(bsu_node) #does not work
print(tree2)

help(tree.prune) #prunes a terminal calde from the tree
# prune is in Bio.Phylo.BaseTree:
# to prune a node using label, I probably can find this terminal node first and then call prune()



# Example Newick and Nexus files
EX_NEWICK = 'Nexus/int_node_labels.nwk'
EX_NEXUS = 'Nexus/test_Nexus_input.nex'

# Example PhyloXML files
EX_APAF = 'PhyloXML/apaf.xml'
EX_BCL2 = 'PhyloXML/bcl_2.xml'
EX_PHYLO = 'PhyloXML/phyloxml_examples.xml'

 





