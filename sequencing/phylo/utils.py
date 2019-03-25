import os
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


def to_hex(n):
    if n == 0:
        return '00'
    return hex(int(n*255))[2:].upper()


def hex_to_rgb(color_map):
    rgb = '#'
    for i in color_map:
        rgb += to_hex(i)
    return rgb


from datetime import datetime
def fix_directories(base_dir, ind_name):
    """
    make sure the dir base_path/curr_date/ind_name exists and create otherwise
    """
    curr_date = "{:%d_%m_%Y}".format(datetime.now())
    dirs_to_create = "{}/{}/{}".format(base_dir,curr_date,ind_name)
    if not os.path.exists(dirs_to_create):
        os.makedirs(dirs_to_create)
    return dirs_to_create


def get_cells_group_map(srs):
    cells_group_map = dict()
    for sr in srs:
        ac = sr.barcoded_content.subclass.amplified_content
        cell = ac.cell
        cells_group_map[sr.id] = cell.custom_group_labeling
    return cells_group_map


cell_labeler = lambda sr: 'ID{}'.format(sr.id)


def textualize_d(objective_d):
    d = dict()
    for sr in objective_d:
        if sr == 'root':
            sr_label = 'root'
        else:
            sr_label = cell_labeler(sr)
        for ms in objective_d[sr]:
            d.setdefault(sr_label, dict())['{}_{}'.format(
                ms.repeat_unit_type,
                ms.id
            )] = objective_d[sr][ms]
    return d