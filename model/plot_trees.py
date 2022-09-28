#Plot BoostMEC trees
#This code is modified from the lightgbm source code for plotting/create_tree_digraph()
#https://lightgbm.readthedocs.io/en/latest/pythonapi/lightgbm.create_tree_digraph.html
#https://lightgbm.readthedocs.io/en/latest/_modules/lightgbm/plotting.html

import os
import matplotlib
import pydot
import lightgbm as lgb
from lightgbm.basic import Booster
from lightgbm.compat import (MATPLOTLIB_INSTALLED, GRAPHVIZ_INSTALLED, LGBMDeprecationWarning,
                     range_, zip_, string_type)
from lightgbm.sklearn import LGBMModel
from pathlib import Path
import pandas as pd
import argparse

#Get desired tree index
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--index', type=int, help='BoostMEC tree index, tells plot_trees.py which tree to plot')
args = parser.parse_args()

tree_index = args.index

#Get file directory
script_path = os.path.realpath(__file__)
script_dir = str(Path(script_path).parent)

lgb_mod = lgb.Booster(model_file = script_dir + '/model.txt')

def _check_not_tuple_of_2_elements(obj, obj_name='obj'):
    """Check object is not tuple or does not have 2 elements."""
    if not isinstance(obj, tuple) or len(obj) != 2:
        raise TypeError('%s must be a tuple of 2 elements.' % obj_name)


def _float2str(value, precision=None):
    return ("{0:.{1}f}".format(value, precision)
            if precision is not None and not isinstance(value, string_type)
            else str(value))

mono_dictionary = {'1':'A', '2':'C', '3':'G', '4':'T'}
di_dictionary = {'1':'AA','2':'AC','3':'AG','4':'AT',
                 '5':'CA','6':'CC','7':'CG','8':'CT',
                 '9':'GA','10':'GC','11':'GG','12':'GT',
                 '13':'TA','14':'TC','15':'TG','16':'TT'}

di_25_dictionary = {'1':'AG', '2':'CG', '3':'GG', '4':'TG'}
di_27_dictionary = {'1':'GA', '2':'GC', '3':'GG', '4':'GT'}

new_names = pd.read_csv(script_dir + "/nucleotide_nice_names.csv")
name_dict = dict(new_names.values)

#modified verison of LightGBM _to_graphviz()
def _to_graphviz2(tree_info, show_info, feature_names, precision=3, **kwargs):
    """Convert specified tree to graphviz instance.

    See:
      - https://graphviz.readthedocs.io/en/stable/api.html#digraph
    """
    if GRAPHVIZ_INSTALLED:
        from graphviz import Digraph
    else:
        raise ImportError('You must install graphviz to plot tree.')
    
    def add(root, parent=None, decision=None):
        """Recursively add node or edge."""
        if 'split_index' in root:  # non-leaf
            name = 'split{0}'.format(root['split_index'])
            l_dec = 'yes'
            r_dec = 'no'
            if root['decision_type'] == '<=':
                operator = "&#8804;"
            elif root['decision_type'] == '==':
                operator = "="
            else:
                raise ValueError('Invalid decision type in tree model.')
            
            if feature_names is not None:
                if feature_names[root['split_feature']] in list(name_dict.keys()):
                    label = f"<B>{name_dict[feature_names[root['split_feature']]]}</B> {operator}"
                else:
                    label = f"<B>{feature_names[root['split_feature']]}</B> {operator}"
            else:
                label = 'split_feature_index: {0} '.format(root['split_feature']) + operator
            
            if feature_names[root['split_feature']] == 'X25di':
                label += f"<B>{'||'.join([di_25_dictionary[key] for key in root['threshold'].split('||')])}</B>"
            elif feature_names[root['split_feature']] == 'X27di':
                label += f"<B>{'||'.join([di_27_dictionary[key] for key in root['threshold'].split('||')])}</B>"
            elif 'di' in feature_names[root['split_feature']]:
                label += f"<B>{'||'.join([di_dictionary[key] for key in root['threshold'].split('||')])}</B>"
            elif 'mono' in feature_names[root['split_feature']]:
                label += f"<B>{'||'.join([mono_dictionary[key] for key in root['threshold'].split('||')])}</B>"
            else:
                label += f"<B>{_float2str(root['threshold'], precision)}</B>"
            for info in show_info:
                if info in {'split_gain', 'internal_value'}:
                    label += r'<br/>{0}: {1}'.format(info, _float2str(root[info], precision))
                elif info == 'internal_count':
                    label += r'<br/>{0}: {1}'.format(info, root[info])
            label = f"<{label}>"
            graph.node(name, label=label, shape = "rectangle")
            add(root['left_child'], name, l_dec)
            add(root['right_child'], name, r_dec)
        else:  # leaf
            name = 'leaf{0}'.format(root['leaf_index'])
            label = 'leaf {0}: '.format(root['leaf_index'])
            label += r'\n{0}'.format(_float2str(root['leaf_value'], precision))
            if 'leaf_count' in show_info:
                label += r'\nleaf_count: {0}'.format(root['leaf_count'])
            graph.node(name, label=label, shape = "ellipse")
        if parent is not None:
            graph.edge(parent, name, decision)
    
    graph = Digraph(**kwargs)
    graph.attr("graph", nodesep="0.05", ranksep="0.3", rankdir="LR")
    add(tree_info['tree_structure'])
    
    return graph


def create_tree_digraph2(booster, tree_index=0, show_info=None, precision=3,
                        old_name=None, old_comment=None, old_filename=None, old_directory=None,
                        old_format=None, old_engine=None, old_encoding=None, old_graph_attr=None,
                        old_node_attr=None, old_edge_attr=None, old_body=None, old_strict=False, **kwargs):
    """Create a digraph representation of specified tree.

    Note
    ----
    For more information please visit
    https://graphviz.readthedocs.io/en/stable/api.html#digraph.

    Parameters
    ----------
    booster : Booster or LGBMModel
        Booster or LGBMModel instance to be converted.
    tree_index : int, optional (default=0)
        The index of a target tree to convert.
    show_info : list of strings or None, optional (default=None)
        What information should be shown in nodes.
        Possible values of list items: 'split_gain', 'internal_value', 'internal_count', 'leaf_count'.
    precision : int or None, optional (default=None)
        Used to restrict the display of floating point values to a certain precision.
    **kwargs
        Other parameters passed to ``Digraph`` constructor.
        Check https://graphviz.readthedocs.io/en/stable/api.html#digraph for the full list of supported parameters.

    Returns
    -------
    graph : graphviz.Digraph
        The digraph representation of specified tree.
    """
    if isinstance(booster, LGBMModel):
        booster = booster.booster_
    elif not isinstance(booster, Booster):
        raise TypeError('booster must be Booster or LGBMModel.')
    
    for param_name in ['old_name', 'old_comment', 'old_filename', 'old_directory',
                       'old_format', 'old_engine', 'old_encoding', 'old_graph_attr',
                       'old_node_attr', 'old_edge_attr', 'old_body']:
        param = locals().get(param_name)
        if param is not None:
            warnings.warn('{0} parameter is deprecated and will be removed in 2.4 version.\n'
                          'Please use **kwargs to pass {1} parameter.'.format(param_name, param_name[4:]),
                          LGBMDeprecationWarning)
            if param_name[4:] not in kwargs:
                kwargs[param_name[4:]] = param
    if locals().get('strict'):
        warnings.warn('old_strict parameter is deprecated and will be removed in 2.4 version.\n'
                      'Please use **kwargs to pass strict parameter.',
                      LGBMDeprecationWarning)
        if 'strict' not in kwargs:
            kwargs['strict'] = True
    
    model = booster.dump_model()
    tree_infos = model['tree_info']
    if 'feature_names' in model:
        feature_names = model['feature_names']
    else:
        feature_names = None
    
    if tree_index < len(tree_infos):
        tree_info = tree_infos[tree_index]
    else:
        raise IndexError('tree_index is out of range.')
    
    if show_info is None:
        show_info = []
    
    graph = _to_graphviz2(tree_info, show_info, feature_names, precision, **kwargs)
    
    return graph


lgb._to_graphviz2 = _to_graphviz2
lgb.create_tree_digraph2 = create_tree_digraph2

plot = lgb.create_tree_digraph2(lgb_mod, tree_index = tree_index, show_info = ['internal_value'])

with open('BoostMEC_tree_' + str(tree_index) + '.svg', 'w') as f:
    f.write(plot._repr_image_svg_xml())
