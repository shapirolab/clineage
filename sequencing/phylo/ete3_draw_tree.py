__author__ = 'veronika'

import argparse
import ete3
import math
import os
import sys

from .utils import get_cluster_colors, get_branch_width, get_leaf_order, get_cells_labels, get_bootsrtap_size, \
                   get_cells_colors, get_dup_labels, get_legend


def size_clustering(t, styles, clustering_sizes_file):
    # add width to branches
    branch_width = get_branch_width(clustering_sizes_file)
    for name, groups in branch_width.items():
        nodes = t.search_nodes(name=name)
        assert len(nodes) == 1, nodes
        node = nodes[0]
        for group, pvalue in groups.items():
            width = min(10, -round(math.log10(pvalue)))
            children = node.get_children()
            for child in children:
                if group not in styles[child.name]:
                    styles[child.name][group] = dict()
                styles[child.name][group]["hz_line_width"] = width
            if group not in styles[name]:
                styles[name][group] = dict()
            styles[name][group]["vt_line_width"] = width
            styles[name]['style']["vt_line_width"] = width+2
    return t, styles


def color_clustering(t, ts, styles, clustering_colors_file):
    # add colors to branches
    cluster_colors = get_cluster_colors(clustering_colors_file)
    for name, groups in cluster_colors.items():
        nodes = t.search_nodes(name=name)
        assert len(nodes) == 1, nodes
        node = nodes[0]
        for group, color in groups.items():
            children = node.get_children()
            for child in children:
                if group not in styles[child.name]:
                    styles[child.name][group] = dict()
                styles[child.name][group]["hz_line_color"] = color

                if 'vt_line_width' in styles[child.name][group]:
                    line_width = child.dist*ts.scale-styles[child.name][group]["vt_line_width"]-0.5
                else:
                    line_width = child.dist*ts.scale
                if 'vt_line_width' in styles[node.name][group]:
                    line_width = line_width + styles[node.name][group]["vt_line_width"]

                rf = ete3.faces.RectFace(height=styles[child.name][group]["hz_line_width"],
                                         width=line_width,
                                         fgcolor=color,
                                         bgcolor=color)
                child.add_face(rf, 0, position='float')
                styles[child.name]['style']["hz_line_width"] = 0
            styles[name][group]["hz_line_color"] = color
            if group not in styles[name]:
                styles[name][group] = dict()
            styles[name][group]["vt_line_color"] = color
            styles[name]['style']["vt_line_color"] = color
    return t, ts, styles


def node_check(name, t):
    nodes = t.search_nodes(name=name)
    assert len(nodes) >= 1
    if len(nodes) == 0:
        Warning("The json file contains too many nodes")
        return None
    node = nodes[0]
    return node


def str2bool(v):
    """
    susendberg's function
    :param v:
    :return:
    """
    return v.lower() in ("yes", "true", "t", "1")


def tree_draw(tree_file,
              tree_name=None,
              order_vector_file=None,
              cell_colors_file=None,
              clustering_colors_file=None,
              clustering_sizes_file=None,
              intermediate_node_sizes_file=None,
              intermediate_node_labels_file=None,
              leaf_labels_file=None,
              legend_file=None,
              duplicate_file=None,
              tree_scale='linear',
              tree_rotation=True,
              font_size=7,
              font_legend=7,
              node_size=3,
              scale_rate=None,
              distance_factor=1,
              y_scale=True
              ):

    t = ete3.Tree(newick=tree_file, format=1)
    ts = ete3.TreeStyle()
    if tree_rotation:
        ts.rotation = 90
    ts.show_leaf_name = True
    ts.show_scale = False
    ts.scale = 1
    if tree_name:
        ts.title.add_face(ete3.TextFace(tree_name, fsize=20), column=0)

    styles = {}
    max_dist = 0

    # initialize all nodes and branches
    for n in t.traverse():
        styles[n.name] = dict()
        styles[n.name]['style'] = ete3.NodeStyle()
        styles[n.name]['style']['fgcolor'] = 'black'
        styles[n.name]['style']["vt_line_width"] = 2
        styles[n.name]['style']["hz_line_width"] = 1
        max_dist = max(max_dist, n.dist)
        # print (max_dist)

    # calculate the scale for the tree (log, linear and right size)
    if tree_scale == 'log':
        max_dist = 0

    root = t.get_tree_root()
    last_leaf = root.get_farthest_leaf()
    ts.y_axis['scale_min_value'] = root.dist
    ts.y_axis['scale_max_value'] = last_leaf[1]

    for n in t.traverse():
        if tree_scale == 'log':
            if n == root:
                styles[n.name]['dist'] = 0
            else:
                father_path = 0
                for ancestor in n.get_ancestors():
                    father_path += styles[ancestor.name]['dist']

                dist = math.log10(n.get_distance(root)*distance_factor+1)-father_path
                if dist < 0:
                    dist = 0
                styles[n.name]['dist'] = dist
                max_dist = max(max_dist, dist)

        elif tree_scale == 'linear':
            # if max_dist > 1:
            #     styles[n.name]['dist'] = round(n.dist/max_dist)
            # else:
            styles[n.name]['dist'] = n.dist

    # leaf styles and update distance
    if not scale_rate:
        # scale_rate = max(10, round(1/max_dist))
        scale_rate = ts.scale

    for n in t.traverse():
        if 'dist' in styles[n.name]:
            n.dist = styles[n.name]['dist']*scale_rate
        if not n.is_leaf():
            styles[n.name]['style']["size"] = 0
        else:
            styles[n.name]['style']["size"] = node_size

    # add bootstrap values to the branches (size of the node)
    if intermediate_node_sizes_file and os.path.isfile(intermediate_node_sizes_file):
        bootsrtap_sizes = get_bootsrtap_size(intermediate_node_sizes_file)
        for branch, size in bootsrtap_sizes.items():
            styles[branch]['style']["size"] = size
            styles[branch]['style']['fgcolor'] = 'black'

    # add colors to the leafs
    if cell_colors_file and os.path.isfile(cell_colors_file):
        cells_colors = get_cells_colors(cell_colors_file)
        for name, color in cells_colors.items():
            styles[name]['style']['fgcolor'] = color

    # reorder the tree by pre-proses if possible
    if order_vector_file and os.path.isfile(order_vector_file):
        leaf_order = get_leaf_order(order_vector_file)
        for n in t.traverse('postorder'):
            if n.get_descendants():
                a = ''
                for leaf in n.get_descendants(strategy='postorder'):
                    if leaf.is_leaf():
                        if not a:
                            a = leaf
                b = n.get_descendants(strategy='preorder')[-1]

                if a.is_leaf() and b.is_leaf():
                    if leaf_order[a.name] > leaf_order[b.name]:
                        left, right = n.children
                        n.children = [right, left]

    # add width to branches
    if clustering_sizes_file and os.path.isfile(clustering_sizes_file):
        t, styles = size_clustering(t, styles, clustering_sizes_file)

    # add colors to branches
    if clustering_colors_file and os.path.isfile(clustering_colors_file):
        t, ts, styles = color_clustering(t, ts, styles, clustering_colors_file)

    # add new leaf labels
    if leaf_labels_file and os.path.isfile(leaf_labels_file):
        cells_labels = get_cells_labels(leaf_labels_file)
        ts.show_leaf_name = False
        for name, label in cells_labels.items():
            nodes = t.search_nodes(name=name)
            assert len(nodes) == 1, nodes
            node = nodes[0]
            if name in cells_colors:
                name_face = ete3.faces.TextFace(cells_labels[name], fsize=font_size, fgcolor=cells_colors[name])
            else:
                name_face = ete3.faces.TextFace(cells_labels[name], fsize=font_size)

            name_face.margin_left = 3
            node.add_face(name_face, 0, "aligned")

    # add duplicate tags to nodes
    if duplicate_file and os.path.isfile(duplicate_file):
        dup_labels = get_dup_labels(duplicate_file)
        for name, color in dup_labels.items():
            node = node_check(name, t)
            if not node:
                continue
            dup_face = ete3.faces.TextFace('*', fsize=10, fgcolor=color)
            dup_face.margin_left = 5
            node.add_face(dup_face, column=1)

    # add y-scale to the picture
    if y_scale:
        ts.y_axis['show'] = True
        ts.y_axis['scale_type'] = tree_scale
        ts.y_axis['scale_length'] = int(last_leaf[1] - root.dist + 10)
        print('ts.y_axis[\'scale_length\']: ', ts.y_axis['scale_length'])

    # add legend to the tree
    if legend_file and os.path.isfile(legend_file):
        legend = get_legend(legend_file)
        for mark in list(legend.keys()):
            ts.legend.add_face(ete3.faces.CircleFace(5, legend[mark]), column=0)
            legend_txt = ete3.faces.TextFace(mark, fsize=font_legend)
            legend_txt.margin_left = 5
            ts.legend.add_face(legend_txt, column=1)
        ts.legend_position = 4

    # set all the styles
    for n in t.traverse():
        n.set_style(styles[n.name]['style'])

    return t, ts
