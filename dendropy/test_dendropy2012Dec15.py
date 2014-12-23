import dendropy

dir(dendropy)
help(dendropy.Tree)

tree_str = "[&R] ((A, (B, (C, (D, E)))),(F, (G, H)));"

tree = dendropy.Tree.get_from_string(
        tree_str,
        "newick")
dir(tree)

print("Before:")
print(tree.as_string('newick'))
print(tree.as_ascii_plot())
tree.prune_taxa_with_labels(["A", "C", "G"])
print("After:")
print(tree.as_string('newick'))
print(tree.as_ascii_plot())


#this seems total length
#trees = dendropy.TreeList.get_from_path("pythonidae.random.bd0301.tre", "nexus")
tree_lengths = [tree.length() for tree in trees]
tree_lengths.sort()
crit_index_95 = int(0.95 * len(tree_lengths))
crit_index_99 = int(0.99 * len(tree_lengths))

print("95%% critical value: %s" % tree_lengths[crit_index_95])
print("99%% critical value: %s" % tree_lengths[crit_index_99])



import dendropy

#mle = dendropy.Tree.get_from_path('examples/pythonidae.mle.nex', 'nexus')
mle = dendropy.Tree.get_from_path('qincodeml/tree.nwk', 'newick')
#denfropy does parse multiple labels output by codeml
# this is annoying

treestring = "(Bha #0.0626 : 0.218136, (Bpu#2 #0.0744 : 0.159123, (Bli#2 #0.0744 : 0.072925, (Bam#2 #0.0744 : 0.096945, Bsu#2 #0.0744 : 0.055166)#2 #0.0744 : 0.074617)#2 #0.0744 : 0.045819)#2 #0.0744 : 0.128114, (Bwe#1 #0.0380 : 0.007139, (Bce#1 #0.0380 : 0.007184, (Ban#3 #0.0001 : 0.000003, Bth#1 #0.0380 : 0.000251)#1 #0.0380 : 0.009027)#1 #0.0380 : 0.004363)#1 #0.0380 : 0.106680);"
tree = dendropy.Tree.geet_from_string(streestring, 'newick)

import re



mle_len = mle.length()
for edge in mle.postorder_edge_iter():
    #edge.length = None
    print(edge)

leaf_edges = mle.leaf_edge_iter()
edge = leaf_edges.next()
edge.is_internal() #return False


print(mle.as_string("newick"))


#branch lengths
trees = dendropy.Tree.get_from_path("int_node_labels.nwk", "newick")
