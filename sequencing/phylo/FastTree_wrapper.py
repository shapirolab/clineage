from plumbum import local
import dendropy
import pandas as pd
from sequencing.phylo.utils import fix_directories, add_root_to_dict

symbols = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
symbol_map = {i+7: s for i, s in enumerate(symbols)}
fasttree = local["/home/dcsoft/s/ron/usefull_apps/fasttree_new/FastTree"]
fasttree_trans_noa = fasttree['-trans', '/home/dcsoft/s/Ofir/noa_for_fasttree_scientific.mat']
fasttree_pseudo = fasttree['-pseudo']
fasttree_pseudo_trans_noa = fasttree_pseudo['-trans', '/home/dcsoft/s/Ofir/noa_for_fasttree_scientific.mat']
fasttree_pseudo_trans_unified = fasttree_pseudo['-trans', '/home/dcsoft/s/Ofir/unified_for_fastree_asaf.mat']
fasttree_trans_unifed = fasttree['-trans', '/home/dcsoft/s/Ofir/unified_for_fastree_asaf.mat']


def create_FastTree_fasta_file(td, rldr, ind_name, title, min_group_size, rldr_group):
    td_with_root = add_root_to_dict(td, rldr)
    df = pd.DataFrame.from_dict(td_with_root)
    vs = list(v for vs in td_with_root.values() for v in vs.values())
    for i in range(min(vs), max(vs) + 4):
        if i not in symbol_map:
            symbol_map[i] = '-'

    df.replace(symbol_map, inplace=True)
    df.fillna("-", inplace=True)
    normalized_fasta_dir = fix_directories("/home/dcsoft/s/trees", ind_name)
    normalized_fasta_path = '{}/{}_{}_FastTree_g{}_root_{}.fasta'.format(normalized_fasta_dir, title, ind_name,
                                                                   min_group_size, rldr_group)
    with open(normalized_fasta_path, 'w') as f:
        for sr, row in df.iterrows():
            f.write('>{}\r\n{}\r\n'.format(sr, ''.join(row.values)))

    return normalized_fasta_path


def run_fasttree_binary(newick_output_file, mutations_fasta_path, variation=fasttree_pseudo_trans_noa):
    variation('-out', newick_output_file, mutations_fasta_path)
    good_tree = dendropy.Tree.get_from_path(
        newick_output_file,
        "newick")
    return good_tree


def run_fastree(fasta_path, ind_name, title, min_group_size, rldr_group):
    newick_output_files = []
    normalized_newick_dir = fix_directories("/home/dcsoft/s/trees", ind_name)
    for name, method, kwargs in [
        ('FTPU', 'fasttree_pseudo_trans_unified', {'variation': fasttree_pseudo_trans_unified}),
    ]:
        newick_output_file = '{}/{}_{}_g{}_{}_root_{}.newick'.format(normalized_newick_dir, title, ind_name,
                                                                     min_group_size, method, rldr_group)
        ret_tree = run_fasttree_binary(newick_output_file, fasta_path, list(kwargs.values())[0])
        newick_output_files.append(newick_output_file)

    return newick_output_files

def unroot_FastTree_tree(newick_output_file):
    rooted_tree = dendropy.Tree.get_from_path(newick_output_file, "newick")
    rooted_tree.reroot_at_node(rooted_tree.find_node_with_taxon_label('root'), update_bipartitions=False)
    tree = dendropy.Tree.get_from_string(rooted_tree.as_string('newick', suppress_rooting=False).replace('root', ''), 'newick')
    tree.ladderize(ascending=False)
    return tree

