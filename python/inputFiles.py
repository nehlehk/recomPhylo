import sys
import re
from dendropy import Tree, DnaCharacterMatrix
import dendropy


def test(tree_path , output=None):
    tree = Tree.get_from_path(tree_path, 'newick')
    print("Original tree")
    print(tree.as_string(schema='newick'))
    print(tree.as_ascii_plot(show_internal_node_labels=True))


if __name__ == '__main__':
    tree_path = input('Enter tree path:')
    test(tree_path)