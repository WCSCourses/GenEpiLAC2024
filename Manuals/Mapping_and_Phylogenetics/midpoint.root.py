import sys

from ete3 import Tree

intre = sys.argv[1]

tre = Tree(intre, quoted_node_names=True)

# Calculate the midpoint node
midpoint = tre.get_midpoint_outgroup()

# Set midpoint as outgroup.
tre.set_outgroup(midpoint)

print(tre.write())