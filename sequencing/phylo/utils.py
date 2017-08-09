import dendropy


def normalize_newick(newick_path, normalized_newick_path):
    """
    Read neick tree from file, delete edge weights and write to another file
    Args:
        newick_path: 
        normalized_newick_path: 

    Returns:

    """
    tns = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get_from_path(
            newick_path,
            "newick",
            taxon_namespace=tns)
    for e in tree.edges():
        e.length = None
    tree.write_to_path(normalized_newick_path, 'newick', suppress_rooting=True)
