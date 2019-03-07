#!/usr/bin/env python3

"""Parse results of codeml (.mlc files) to save dN, dS values in a table"""


from __future__ import print_function

import sys
import os.path
import re
import numpy as np
import ete3
import argparse
import logging

import LibsDyogen.myPhylTree as PhylTree # my custom python3 version

from genomicustools.identify import convert_gene2species
from IOtools import Stream

logger = logging.getLogger(__name__)


#from codeml.codemlparser import mlc_parser

np.set_printoptions(formatter={"float_kind": lambda x: "%g" %x})

ENSEMBL_VERSION = 85
PHYLTREEFILE = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data{0:d}/PhylTree.Ensembl.{0:d}.conf"
ANCGENE2SP = re.compile(r'([A-Z][A-Za-z0-9_.-]+)ENS')
# ~~> genomicus.my_identify?


def printtree(tree, indent='', features=None, **kwargs):
    line = indent + tree.name
    if features:
        line += ': ' + ' '.join('%s=%s' % (ft,getattr(tree, ft))
                                for ft in features
                                if getattr(tree, ft, None))
    print(line, **kwargs)
    for ch in tree.children:
        printtree(ch, (indent+'  '), features)


def showtree(fulltree, ages):
    """Default `showtree` function: do nothing. This function can be
    redefined using `def_showtree`"""
    pass



def def_showtree(measures, show=None):  # ~~> genomicus.reconciled
#def def_showreconciledgenetree
    """Depending on the boolean argument 'show', return the function
    `showtree`:
    if show=None, this function does nothing;
    else, showtree display an interactive tree using ete3."""
    if show:
        if show == 'gui':
            def show_func(tree, ts=None):
                print("Showing tree.")
                tree.show(tree_style=ts, name=tree.name)
        elif show == 'notebook':
            def show_func(tree, ts=None):
                #print("Rendering tree.", file=sys.stderr)
                print(("Defining global FULLTREE and TS. Display with:\n"
                       "    FULLTREE.render('%%inline', tree_style=TS)"))
                #tree.render('%%inline', w=500, tree_style=ts) # 
                global FULLTREE, TS
                FULLTREE = tree
                TS = ts
        else:
            # Save to an image file.
            def show_func(tree, ts=None):
                print("Saving tree in %r." % show)
                tree.render(show, tree_style=ts)

        # define duplication node style:
        ns = {ntype: ete3.NodeStyle() for ntype in ('spe', 'dup', 'leaf', 'root')}
        ns['dup']['fgcolor'] = 'red'
        ns['leaf']['fgcolor'] = 'green'
        ns['root']['fgcolor'] = 'black'
            
        header = measures + ['age_'+m for m in measures] + ['type']
        if 'dist' in header: header.remove('dist')
        measure_faces = {feature: ete3.AttrFace(feature,
                                                text_prefix='%s: ' % feature,
                                                fsize=7,
                                                fgcolor='grey',
                                                formatter='%.4g')
                            for feature in header[:-1]}

        def mylayout(node):
            for feature in header[:-1]:
                if hasattr(node, feature):
                    node.add_face(measure_faces[feature], column=0,
                                  position='branch-bottom')
            node.set_style(ns[getattr(node, 'type', 'root')])

        #leaf_info = dict(zip(header, [np.NaN]))

        def showtree(fulltree):
            # TODO: this should be a layout function. Then change ts.layout_fn.
            # define a tree style:
            ts = ete3.TreeStyle()
            ts.show_branch_length = True
            ts.layout_fn = mylayout

            show_func(fulltree, ts)
            return fulltree
    else:
        def showtree(fulltree):
            pass
    return showtree


def load_fulltree(mlcfile, replace_nwk='.mlc', replace_by='.nwk'):
    """Return the ete3.Tree object corresponding to a given .mlc file.

    It simply finds the newick file by replacing the extension of the mlc file.
    Catch errors (file does not exist / wrong format)"""
    #nwkfile = mlcfile.replace(replace_nwk, '.nwk')
    rootname = re.sub(re.escape(replace_nwk) + '$', '', mlcfile)
    nwkfile = rootname + replace_by
    try:
        fulltree = ete3.Tree(nwkfile, format=1)
    except ete3.parser.newick.NewickError as e:
        if os.path.exists(nwkfile):
            logger.error("\nNewickError: Malformed newick tree structure in %r",
                          nwkfile)
        else:
            logger.error("\nNewickError: Unexisting tree file %r" % nwkfile)
        sys.exit(2)

    fulltree.add_feature('treename', os.path.basename(rootname))
    return fulltree


def eat_lines_uptomatch(file_obj, regex, line=None, jump=0):
    if line is None:
        line = file_obj.readline()
    match = regex.match(line.rstrip())
    while not match:
        for _ in range(jump):
            file_obj.readline()

        line = file_obj.readline()
        if not line:
            raise LookupError("Reached end of file without matching '%s'" % regex.pattern)
        match = regex.match(line.rstrip())
    return line #, match


def yield_lines_whilematch(file_obj, regex, line=None, jump=0):
    if line is None:
        line = file_obj.readline()
    match = regex.match(line.rstrip())
    while match:
        yield line
        for _ in range(jump):
            file_obj.readline()

        line = file_obj.readline()
        match = regex.match(line.rstrip())


def branch2nb(mlc, fulltree):  # ~~> codeml.codeml_parser?
    """Parse the codeml result file (.mlc) to return 2 trees:
    tree_nbs: the tree with node labelled as numbers.
    tree_ids: the tree with original node labels (including inner nodes).
    
    Also return the dictionaries to convert nb -> id and conversely.

    Arguments:
        - mlc     : an opened file (codeml result file)
        - fulltree: the ete3.Tree for the gene tree (with reinserted nodes)
    """
    regex_tree = re.compile(r'^(.*);$')
    regex_lnL  = re.compile(r'^lnL\(')
    regex_w    = re.compile(r'^w \(dN/dS\) for branches:')
    regex_freqtable = re.compile('^Codon position x base \(3x4\) table for each sequence\.$')
    regex_seqid = re.compile('^#(\d+): (.*)$')

    line = eat_lines_uptomatch(mlc, regex_freqtable)

    seqids = {}

    mlc.readline()
    line = mlc.readline()
    while True:
        m = regex_seqid.match(line.rstrip())
        if not m:
            break
        seqids[m.group(1)] = m.group(2)
        for _ in range(6):
            line = mlc.readline()
        assert line
        #if not line:
        #    break

    assert seqids
    
    # Get the line listing all branches
    while not regex_lnL.match(line):
        line = mlc.readline().rstrip()
    branches_line = mlc.readline().rstrip()
    branches = branches_line.split()
    # get translation by looking at the newick tree lines.
    lengths_line = mlc.readline().rstrip()
    lengths = [float(l) for l in lengths_line.split()]

    # Get the tree labelled with numbers
    line = mlc.readline().rstrip()
    while not regex_tree.match(line):
        line = mlc.readline().rstrip()
    tree_nbs = ete3.Tree(line)

    # Get the tree with original labels
    line = mlc.readline().rstrip()
    while not regex_tree.match(line):
        line = mlc.readline().rstrip()
    tree_ids = ete3.Tree(line)

    # Get the omega (dN/dS) values
    line = mlc.readline().rstrip()
    while not regex_w.match(line):
        line = mlc.readline().rstrip()
    omegas = [float(w) for w in regex_w.sub('', line).split()]
    branch_tw = list(zip(branches, lengths, omegas))

    id2nb = dict(zip(tree_ids.get_leaf_names(), tree_nbs.get_leaf_names()))
    nb2id = dict(zip(tree_nbs.get_leaf_names(), tree_ids.get_leaf_names()))

    assert set(seqids) == set(nb2id) and all((seqid == seqids[nb]) for nb,seqid in nb2id.items())

    # Remember nodes that were detached from fulltree, and all their descendants
    detached_subtrees = set()

    # get internal nodes nb-to-id conversion (fulltree: tree with internal node
    # annotations)
    for leafid in fulltree.get_leaves():
        try:
            leafid.add_feature('nb', id2nb[leafid.name])
        except KeyError as err:
            # This leaf is not a species node. skip.
            pass

    # branches follow the order of the newick string.
    while branches:
        br = branches.pop()
        base, tip = br.split('..')
        debug_msg = "%-8s " % br
        # Update the tree_nbs
        base_nb_node = tree_nbs.search_nodes(name=tip)[0].up
        base_nb_node.name = base
        # Then find the matching node in fulltree
        try:
            tip_id = nb2id[tip]
            found_tips = fulltree.search_nodes(name=tip_id)
            try:
                base_node = found_tips[0].up
            except IndexError as err:
                logger.warning('Node %s:%r not found in fulltree', tip, tip_id)
                      #file=sys.stderr)
                # TODO: search in tree_ids, then detached_subtrees.add()
                detached_subtrees.add(tip)
                continue

                #while len(base_node.children) == 1:
                #    base_node = base_node.up
            except AttributeError as err:
                logger.warning('root (%r: %r) cannot have ancestors', tip, tip_id)
                #print(fulltree.search_nodes(name=tip_id), err)
                continue

            base_id = base_node.name
            debug_msg += "%s -> %s  " % (base_id, tip_id)
            nb2id[base] = base_id
            id2nb[base_id] = base
            # Add number in the fulltree:
            base_node.add_feature('nb', base)
            logger.debug(debug_msg + 'Ok')
        except KeyError as e:
            #if base in detached_subtrees:
            # I assume a progression from leaves to root, otherwise I will miss nodes
            logger.debug(debug_msg + 'Detached')
            #else:
            # Not found now, put it back in the queue for later
            #    branches.insert(0, (base, tip))
            #    logger.debug('KeyError')

    logger.debug(tree_nbs.get_ascii())
    logger.debug(tree_ids.get_ascii())
    logger.debug(fulltree)
    #logger.debug(fulltree.get_ascii())
    #printtree(fulltree, features=['nb'])
    #showtree(fulltree)
    return id2nb, nb2id, tree_nbs, branch_tw

# Needs a unit test for the above (are the missing nodes properly re-inserted?)

def get_dNdS(mlc):  # ~~> codeml.codeml_parser?
    """Parse table of dN/dS from codeml output file.
    
    mlc: filehandle

    WARNING: Numbers are rounded!!! (e.g, only 3 decimals for `t`, 4 for `dS`).
    """
    #reg_dNdS = re.compile(r'dN & dS for each branch')
    reg_dNdS = re.compile(r'\s+branch\s+t\s+N\s+S\s+dN/dS\s+dN\s+dS\s+N\*dN' \
                          r'\s+S\*dS')
    reg_dStree = re.compile(r'dS tree:$')
    reg_dNtree = re.compile(r'dN tree:$')
    
    line = mlc.readline().rstrip()
    while not reg_dNdS.match(line):
        line = mlc.readline().rstrip()
    assert mlc.readline().rstrip() == ''  # skip blank line
    dNdS = {}
    dNdS['colnames'] = line.split()[1:]
    line = mlc.readline().rstrip()
    while line != '':
        #print line
        fields = line.split()
        dNdS[fields[0]] = [float(x) for x in fields[1:]]
        line = mlc.readline().rstrip()

    line = mlc.readline().rstrip()
    while not reg_dStree.match(line):
        line = mlc.readline().rstrip()
    dStree = mlc.readline().rstrip()
    
    line = mlc.readline().rstrip()
    while not reg_dNtree.match(line):
        line = mlc.readline().rstrip()
    dNtree = mlc.readline().rstrip()
    
    return dNdS, dStree, dNtree
    

def tree_nb_annotate(tree, id2nb, tree_nbs):  # ~~> codeml.codeml_parser?
    """Add internal node names (numbers used by codeml) in the tree structure."""
    parent_nbs = {tuple(sorted(ch.name for ch in n.children)): n.name \
                for n in tree_nbs.traverse() if not n.is_leaf()}
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            node.add_feature('nb', id2nb[node.name])
        else:
            node.add_feature('nb',
                    parent_nbs[tuple(sorted(ch.nb for ch in node.children))])


def dNdS_precise(dNdS, br_tw, dStree, dNtree, id2nb, tree_nbs):  # ~~> codeml_parser?
    """Make the dNdS table from the codeml output file more precise:
    dS and dN values have more decimal in the tree string than in the table."""
    # First, update the dNdS table with the more accurate values br_lengths
    # (br_len is actually 't' in the table)
    colindex = {colname: i for i, colname in enumerate(dNdS['colnames'])}
    for br, br_len, br_w in br_tw:
        dNdS[br][colindex['t']] = br_len # 't'
        dNdS[br][colindex['dN/dS']] = br_w   # 'dN/dS'

    # Then update the dNdS table with the more accurate values in dStree dNtree
    dStr = ete3.Tree(dStree)
    dNtr = ete3.Tree(dNtree)
    tree_nb_annotate(dStr, id2nb, tree_nbs)
    tree_nb_annotate(dNtr, id2nb, tree_nbs)
    for col, tr in ((colindex['dS'], dStr), (colindex['dN'], dNtr)):
        for node in tr.traverse():
            if not node.is_leaf():
                for child in node.children:
                    br = '%s..%s' % (node.nb, child.nb)
                    dNdS[br][col] = child.dist
    return dNdS


#def set_dS_fulltree(fulltree, id2nb, dNdS):
def set_dNdS_fulltree(fulltree, id2nb, dNdS, raise_at_intermediates=True):
    """Add codeml dS on each branch of the complete tree (with missing species
    nodes).
    When a node has a single child: take the branch dS on which it is, 
    and divide proportionnally to each segment dist."""

    colindex = {colname: i for i, colname in enumerate(dNdS['colnames'])}

    for node in fulltree.get_descendants('postorder'): # exclude root
        if len(node.children) != 1:
            parent = node.up
            intermediates = []
            #at root, up is None
            while parent.up and getattr(parent, 'reinserted', None):
                intermediates.append(parent)
                parent = parent.up
            #print("%s:%s -> %s:%s\n   %s" % (id2nb[parent.name],
            #                          parent.name,
            #                          id2nb[node.name],
            #                          node.name,
            #                          [ch.name for ch in parent.children]), file=sys.stderr)
            
            if parent.is_root():
                # Disable dating of nodes immediately below the root, as the
                # tree should be considered unrooted and those values aren't
                # reliable.
                logger.info("Node below root: %s. Discard measures.", node.name)
                totals = {valname: np.NaN for valname in colindex}
            else:
                try:
                    br = id2nb[parent.name] + '..' + id2nb[node.name]
                    totals = {valname: dNdS[br][i] for valname, i in colindex.items()}

                except KeyError as err:
                    raise
    
            if intermediates: # this `if` is not necessary, maybe more efficient
                if raise_at_intermediates:
                    raise AssertionError('Node %r has a single child.' % parent.name)

                dist_tot = sum(n.dist for n in (intermediates + [node]))

                if dist_tot == 0:
                    logger.warning("WARNING: dS = 0 between %r and %r",
                                    parent.name, node.name)
                    dist_tot = 1. # should set to NaN to be sure.

                # Constant rates extrapolation.
                for inter_node in (intermediates + [node]):
                    dist_ratio = inter_node.dist / dist_tot
                    for valname, val in totals.items():
                        if valname in ('t', 'dN', 'dS', 'N*dN', 'S*dS'):
                            inter_node.add_feature(valname, val * dist_ratio)
                        else:
                            inter_node.add_feature(valname, val)
            else:
                for valname, val in totals.items():
                    node.add_feature(valname, val)


def rm_erroneous_ancestors(fulltree, phyltree):  # ~~> genomicustools
    """Some ancestors with only one child are present in the species phylogeny.
    They must be removed when they have an age of zero"""
    infinite_dist = 100
    for node in fulltree.iter_descendants():
        if node.dist > infinite_dist:
            logger.warning('DETACH node with dist>%d : %s, dist=%d',
                  infinite_dist, node.name, node.dist)
            # Detaching/deleting would create a mismatch between fulltree and
            # tree_nbs/tree_ids
            #node.add_feature('skip', True)
            parent = node.up
            node.detach()

        if not node.is_leaf():
            try:
                taxon = ANCGENE2SP.match(node.name).group(1).replace('.', ' ')
            except AttributeError:
                raise ValueError("Can not match species name in %r" % \
                                   node.name)
            if len(node.children) <= 1:
                if hasattr(node, 'reinserted'):
                    if phyltree.ages[taxon] > 0:
                        logger.warning('DELETE reinserted node with age>0: %s',
                                        node.name)
                        ### TODO: see if I should actually keep these nodes
                        ###       (most often known ancestors, so might be
                        ###        useful)
                else:
                    logger.warning('DELETE node with single child: %s',
                                    node.name)
                node.delete(prevent_nondicotomic=False,
                            preserve_branch_length=True)


def sum_average_dNdS(dNdS, nb2id, tree_nbs):
    """Compute the average dS depth of each node (between this node and
    leaves).
    
    Probably the roughest method. Use rather `bound_average_dS`"""
    dS, dN = {}, {}
    for n in tree_nbs.traverse('postorder'):
        node = n.name
        if n.is_leaf():
            dS[nb2id[node]] = 0
            dN[nb2id[node]] = 0
            dS[node] = 0
            dN[node] = 0
        else:

            dS_branch = dNdS[node + '..' + c.name][5]
            dN_branch = dNdS[node + '..' + c.name][4]
            dS_children = [dS[nb2id[c.name]] + dS_branch for c in n.children]
            dN_children = [dN[nb2id[c.name]] + dN_branch for c in n.children]
            #dS_children = [dS[c.name]] for c in n.children]
            #dN_children = [dN[c.name]] for c in n.children]
            average_dS_ch = sum(dS_children) / len(dS_children)
            average_dN_ch = sum(dN_children) / len(dN_children)

            dS[nb2id[node]] = average_dS_ch
            dN[nb2id[node]] = average_dN_ch
            dS[node] = average_dS_ch
            dN[node] = average_dN_ch

    return dN, dS


def rec_average_u(node, scname, subtree, measures=['dS'], weightkey='cal_leaves'):
    """Perform the averaging of one node, given its descendants (unweighted).
    
    This operation must be repeated along the tree, from the leaves to the root.

    Update the subtree dictionary (which holds the temporary measures for nodes).
    """
    # weightkey could also be: age

    # Array operation here: increment tmp_m with current measure.
    # One child per row.
    children_ms = np.stack([getattr(ch, m, np.NaN) for m in measures] \
                            + subtree[ch.name]['tmp_m']
                           for ch in node.children)
    # weights:
    children_ws = np.array([[subtree[ch.name][weightkey]] for ch in
                            node.children])
    subtree[scname][weightkey] = children_ws.sum()
    children_ms *= children_ws / children_ws.sum()

    # mean of each measure (across children)
    subtree[scname]['tmp_m'] = children_ms.sum(axis=0)
    #return children_ms.sum(axis=0)


def rec_average_w(node, scname, subtree, measures=['dS']): #, weightkey
    """Perform the averaging of one node, given its descendants (weighted).
    
    This operation must be repeated along the tree, from the leaves to the root
    """
    # Array operation here: increment tmp_m with current measure.
    # One child per row.
    children_ms = np.stack([getattr(ch, m, np.NaN) for m in measures] \
                            + subtree[ch.name]['tmp_m']
                           for ch in node.children)
    # mean of each measure (across children)
    subtree[scname]['tmp_m'] = children_ms.mean(axis=0)
    #return children_ms.mean(axis=0)


# Functions used to determine whether a node should be calibrated in `bound_average`
def isdup(node, taxon, subtree):
    """Checks whether node is a duplication.

    Also see `genetree_drawer.infer_gene_event`.
    """
    children_taxa = set((subtree[ch.name]['taxon'] for ch in 
                         node.children))
    return len( children_taxa & set((taxon,)) ) == 1

def isinternal(node, taxon=None, subtree=None):
    """True if node is neither a leaf nor the root"""
    return not node.is_leaf() and not node.is_root()

#def is_descendant_of(node, taxon=None, subtree=None):

# Example of a check for a given speciation node. Useless because cannot calibrate dups
#def is_murinae_spe(node, taxon, subtree):
#    return taxon=='Murinae' and not isdup(node, taxon, subtree)

def def_is_target_speciation(target_sp):
    """define the speciation node check function"""
    def is_target_spe(node, taxon, subtree):
        return taxon==target_sp and not isdup(node, taxon, subtree)
    return is_target_spe

def def_is_any_taxon(*target_taxa):
    """define the speciation node check function"""
    def is_any_taxon(node, taxon, subtree):
        return taxon in target_taxa
    return is_any_taxon

def negate(testfunc):
    return lambda *args,**kwargs: not testfunc(*args, **kwargs)

def get_tocalibrate(key):
    """'Parse' the key argument to combine several tests together.
    
    Syntax:
        - any whitespace means AND.
        - '|' means OR. AND takes precedence over OR.
        - tests may be any of: 'isdup', 'isinternal', 'taxon'.
           The 'taxon' test must be written with taxa as arguments:
               'taxon:mytaxon1,mytaxon2'.
        - tests can be negated by prefixing with '!'.

    Example to date all internal nodes except a calibrated speciation:
    'isdup | isinternal !taxon:Simiiformes'
    """
    tocalibrate_funcs = {'isdup': isdup, 'd': isdup,
            'isinternal': isinternal, 'isint': isinternal, 'i': isinternal,
            'taxon': def_is_any_taxon, 't': def_is_any_taxon}

    ## TODO: parse parentheses for precedence.
    alt_tests = [] # Split by 'OR' operators first.
    for alt_teststr in key.split('|'):
        tests = []
        for teststr in alt_teststr.split():
            applynegate = False
            if teststr[0] == '!':
                teststr = teststr[1:]
                applynegate = True
            if ':' in teststr:
                teststr, testargs = teststr.split(':')
                testargs = testargs.split(',')
                test = tocalibrate_funcs[teststr](*testargs)
            else:
                test = tocalibrate_funcs[teststr]

            if applynegate:
                test = negate(test)

            tests.append(test)

        alt_tests.append(lambda *args: all(test(*args) for test in tests))

    return lambda *args: any(alt_test(*args) for alt_test in alt_tests)
    

def bound_average(fulltree, ensembl_version, calibration, measures=['dS'],
                    unweighted=False, method2=False,
                    tocalibrate="isdup", keeproot=False):
    """normalize duplication node position between speciation nodes.
     
     Scaling:
     --------
       method1: d / (D' + d) with D' modified (UPGMA/WPGMA)
       method2: d / (D' + d) with D' as the original measure.
        - d: average (WPGMA/UPGMA-like) of branch length leading to next
             speciation, in units of 'measure'.
             'measure' must be an attribute present in fulltree ('dS', 'dist').
        - D': branch length from the previous speciation to this duplication.
          
    Arguments:
    ----------
      - calibration: dictionary of taxon -> age
      - tocalibrate: function returning True when the node should be calibrated.
                     Takes 2 params: (node, taxon).
                     By default: check whether node is a duplication."""

    rec_average = rec_average_u if unweighted else rec_average_w

    # Convert argument to function
    tocalibrate = get_tocalibrate(tocalibrate)

    ages = []
    # list of ete3 subtrees of fulltree (between two consecutive speciations)
    subtrees = []
    subtree = {} # temporary subtree while traversing from one speciation to
    # another
    n_measures = len(measures)
    # initial measure while passing a speciation for the first time
    measures_zeros = np.zeros(n_measures)
    is_subtree_leaf = lambda node: (not node.is_root() and node.cal is True)

    for node in fulltree.traverse('postorder'):
        scname = node.name
        debug_msg = "* %s: " % scname
        ### DISCARD all dup just after the root speciation/duplication.
        if node.is_root() and not keeproot:
            logger.debug(debug_msg + "Root (discard)")
            continue
        #try:
        branch_measures = np.array([getattr(node, m, np.NaN) for m in measures])
        #except AttributeError as err:
        #    if not node.is_root():
        #        err.args += ("Node: %r (is_root: %s)" % (node, node.is_root()),)
        #        raise
        if node.is_leaf():
            eventtype = 'leaf'
            node.add_features(type='leaf', cal=True)
            try:
                taxon = convert_gene2species(scname, ensembl_version)
                leaf_age = 0 # TODO: take From the calibration dictionary
            except KeyError as err:
                logger.warning(err)
                taxon = None
                leaf_age = np.NaN

            subtree[scname] = {'taxon': taxon,
                               'tmp_m': measures_zeros,
                               'br_m': branch_measures,
                               'age': leaf_age,
                               'cal_leaves': 1} # number of descendant calibrated nodes
            #ages[scname] = 0
            ages.append([scname] + branch_measures.tolist() +
                        measures_zeros.tolist() +
                        [1, "leaf", getattr(node.up, 'name', None), taxon,
                            fulltree.name, fulltree.treename])
            logger.debug(debug_msg + "Leaf")
        else:
            #try:
            #    assert len(node.children) > 1
            #except AssertionError as err:
            #    if node.is_root():
            #        pass
            #    else:
            #        raise
            try:
                taxon = ANCGENE2SP.match(scname).group(1).replace('.', ' ')
            except AttributeError:
                raise ValueError("Can not match species name in %r" % scname)

            subtree[scname] = {'taxon': taxon, 'br_m': branch_measures}

            # Compute the temporary measure of the node.
            rec_average(node, scname, subtree, measures)

            eventtype = 'dup' if isdup(node, taxon, subtree) else 'spe'
            node.add_feature('type', eventtype)

            if tocalibrate(node, taxon, subtree):
                # it is uncalibrated:
                node.add_feature('cal', False)
                logger.debug(debug_msg + "Uncal.; m=%s", subtree[scname]['tmp_m'])

                # Store the age of the **next calibrated** node (propagate down).
                # Since it is a duplication, this should be the same for
                # children[0] and children[1]:
                ch_ages = [subtree[ch.name]['age'] for ch in node.children]
                if len(set(ch_ages)) > 1:

                    logger.warning("At %r: unequal children's ages (next "
                                    "calibrated ages): %s", scname, ch_ages)
                    #showtree(fulltree)
                    subtree[scname]['age'] = np.NaN
                else:
                    # This 'age' will *later* be modified (rescaled according to dS)
                    subtree[scname]['age'] = ch_ages[0]

            else:
                # it is calibrated.
                node.add_feature('cal', True)
                logger.debug(debug_msg + "Calibrated")
                # store the age of this taxon
                node_age = subtree[scname]['age'] = calibration[taxon]
                subtree[scname]['cal_leaves'] = 1
                ages.append([scname] + branch_measures.tolist() +
                            [node_age] * n_measures +
                            [1, # calibrated
                             eventtype, getattr(node.up, 'name', None), taxon,
                             fulltree.name, fulltree.treename])

                
                node.add_features(**{'age_'+m: node_age for m in measures})
                
                # Climb up tree until next speciation and assign an age to
                # each duplication met on the way

                # age between this node and the descendant calibrations
                dists_to_calib = [node_age - subtree[ch.name]['age'] \
                                  for ch in node.children]
                mean_dists_to_calib = sum(dists_to_calib) / len(dists_to_calib)

                for ch, dist_to_calib in zip(node.children, dists_to_calib):
                    ch_taxon = subtree[ch.name]['taxon']
                    ### This is where version _2 is different
                    # walk descendants until speciation.
                    # need to get the age of next speciation and compute the
                    # time between the two speciation.

                    if not method2:
                        ### !!! Method1 does not make sense (both sides don't end
                        ###     at the same speciation date).
                        ###     Or at least for each child, choose a proportion of
                        ###     scaling_m corresponding to the relative speciation
                        ###     distance. TODO.
                        scaling_m = subtree[scname]['tmp_m'] * dist_to_calib/mean_dists_to_calib
                    else:
                        scaling_m = None

                    logger.debug("    climb up to next calibration: " \
                                     "dist_to_calib=%s scaling_m=%s tmp_m=%s" \
                                     % (dist_to_calib, scaling_m,
                                        subtree[scname]['tmp_m']))
                    # Dynamic list during the climbing:
                    # store (node, measures for this node's branch)
                    nextnodes = [(ch, subtree[ch.name]['br_m'])]
                    while nextnodes:
                        nextnode, next_path_m = nextnodes.pop(0)
                        try:
                            nextnode_m = subtree[nextnode.name]['tmp_m']
                        except KeyError as err:
                            err.args += ("Error: Node exists twice in the "
                              "tree. You may need to rerun `prune2family.py`",)
                            raise
                        logger.debug("    - %s: measure(to calib)=%s"\
                                         "; measure(from calib)=%s" % \
                                            (nextnode.name,
                                             nextnode_m, next_path_m))
                        # or: any(nextnode_m > measures_zeros) ?
                        if nextnode_m is not measures_zeros:
                            # It's a leaf
                            if method2:
                                scaling_m = next_path_m + nextnode_m

                            if any(scaling_m == 0):
                                logger.warning("Scaling measure = %s ("
                                                "cannot divide) at %r",
                                                scaling_m, nextnode.name)
                            age = node_age - \
                                  (1 - nextnode_m/scaling_m) * dist_to_calib
                            #ages[nextnode.name] = age
                            nextnode_measures = subtree[nextnode.name]['br_m'].tolist()
                            ages.append([nextnode.name] + nextnode_measures + \
                                        age.tolist() + \
                                        [0, #dated
                                         nextnode.type, nextnode.up.name,
                                         subtree[nextnode.name]['taxon'],
                                         #ch_taxon,
                                         fulltree.name,
                                         fulltree.treename])
                            for i, m in enumerate(measures):
                                nextnode.add_feature('age_'+m, age[i])
                            # When climbing up to next spe, need to increment
                            # the lengths 'next_path_m' from the previous spe
                            nextnodes.extend(
                                    (nch, subtree[nch.name]['br_m'] + next_path_m)
                                    for nch in nextnode.children)
                        subtree.pop(nextnode.name)

                # copy the node so that it becomes a rooted tree
                nodecopy = node.copy()
                cal_leaves = nodecopy.get_leaves(is_leaf_fn=is_subtree_leaf)
                for clf in cal_leaves:
                    for clfch in clf.children:
                        clfch.detach()
                
                subtrees.append(nodecopy)

                # then reset measure (dS, dist) to zero
                subtree[scname]['tmp_m'] = measures_zeros
    #logger.debug(subtree)
    return ages, subtrees


def save_ages(ages, opened_outfile):
    for row in ages:
        opened_outfile.write("\t".join(str(elem) for elem in row)
        #                                if not np.isnan(elem) else '')
                             + "\n")

def del_singletons(tree):
    for node in tree.iter_descendants('postorder'):
        #print(node.name, node.children)
        if len(node.children) == 1:
            child = node.children[0]
            #parent = node.up
            for feature in ['branch_dS', 'branch_dN', 'branch_t']:
                if hasattr(child, feature):
                    setattr(child, feature,
                            getattr(child, feature) + getattr(node, feature, np.nan))
            logger.info('DELETE %r before saving', node.name)
            node.delete(prevent_nondicotomic=False, preserve_branch_length=True)
    while len(tree.children) == 1:
        tree = tree.children[0].detach()

    return tree


def save_fulltree(fulltree, opened_outfile):
    # first delete nodes with single child, because R Ape cannot read it.
    #print("Save_fulltree")
    #print(fulltree.get_ascii())
    fulltree = del_singletons(fulltree)  # Move that outside the script (unix philosophy)

    if fulltree.children:
        opened_outfile.write(fulltree.write(features=['cal', 'type',
                                                      'branch_dist',
                                                      't',
                                                      'dS',
                                                      'dN',
                                                      'age_dist',
                                                      'age_t',
                                                      'age_dS',
                                                      'age_dN'],
                                            format=1,
                                            format_root_node=True) + '\n')

def save_subtrees(subtrees, opened_outfile):
    #for stree in subtrees:
    #    opened_outfile.write(stree.write(features=['type'], format=1) + '\n')
    #opened_outfile.write('\n'.join(subtrees) + '\n')
    for stree in subtrees:
        save_fulltree(stree, opened_outfile)

savefunctions = {'ages': save_ages,
                 'fulltree': save_fulltree,
                 'subtrees': save_subtrees} #save_subtrees}


def setup_fulltree(mlcfile, phyltree, replace_nwk='.mlc', replace_by='.nwk',
                   measures=['dS']):
    fulltree = load_fulltree(mlcfile, replace_nwk, replace_by)
    rm_erroneous_ancestors(fulltree, phyltree)
    if set(('dN', 'dS', 't')) & set(measures):
        with open(mlcfile) as mlc:
            id2nb, nb2id, tree_nbs, br_tw = branch2nb(mlc, fulltree)
            dNdS, dStree, dNtree = get_dNdS(mlc)

        dNdS = dNdS_precise(dNdS, br_tw, dStree, dNtree, id2nb, tree_nbs)
        set_dNdS_fulltree(fulltree, id2nb, dNdS)
    return fulltree


#def process_tonewick(mlcfile, phyltree, replace_nwk='.mlc', measures=['dS'], 
#                     **kwargs):
#    fulltree = setup_fulltree(mlcfile, phyltree, replace_nwk, measures)
#
#    # replace the branch distances by the chosen measure
#    measure = measures[0]
#    if measure != 'dist':
#        for node in fulltree.traverse():
#            node.add_feature('old_dist', node.dist)
#            node.dist = getattr(node, measure, np.NaN)
#
#    return fulltree

def set_subtrees_distances(subtrees, measure):
    """Set the new branch lengths to the given measure (dS, dN or t)"""
    for subtree in subtrees:
        for node in subtree.traverse():
            node.add_feature("branch_dist", node.dist)
            newvalue = getattr(node, measure)
            node.dist = 0 if np.isnan(newvalue) else newvalue


#def process_toages(mlcfile, phyltree, replace_nwk='.mlc', measures=['dS'],
def process(mlcfile, ensembl_version, phyltree, replace_nwk='.mlc', replace_by='.nwk',
            measures=['dS'], method2=False, unweighted=False, tocalibrate="isdup", 
            keeproot=False):
    fulltree = setup_fulltree(mlcfile, phyltree, replace_nwk, replace_by, measures)
    
    ages, subtrees = bound_average(fulltree, ensembl_version, phyltree.ages,
                                   measures=measures, unweighted=unweighted,
                                   method2=method2, tocalibrate=tocalibrate, 
                                   keeproot=keeproot)
    showtree(fulltree)
    set_subtrees_distances(subtrees, measures[0])
    return ages, fulltree, subtrees


def main(outfile, mlcfiles, ensembl_version=ENSEMBL_VERSION,
         phyltreefile=PHYLTREEFILE, method2=False,
         measures=['t', 'dN', 'dS', 'dist'], unweighted=False, verbose=False,
         show=None, replace_nwk='.mlc', replace_by='.nwk', ignore_errors=False,
         saveas='ages', tocalibrate='isdup', keeproot=False):
    nb_mlc = len(mlcfiles)
    
    loglevel = logging.DEBUG if verbose else logging.WARNING
    logger.setLevel(loglevel)

    logging.debug("outfile  %s\n"
                  "     mlcfiles      %s %s\n"
                  "     ensembl       v.%d\n"
                  "     measures      %s\n"
                  "     verbose       %s\n"
                  "     show          %s\n"
                  "     method2       %s\n"
                  "     unweighted    %s\n"
                  "     replace_nwk   %s\n"
                  "     replace_by    %s\n"
                  "     ignore_errors %s\n",
                  outfile,
                  (mlcfiles[:5], '...' * (nb_mlc>5)),
                  ensembl_version,
                  measures,
                  verbose,
                  show,
                  method2,
                  unweighted,
                  replace_nwk,
                  replace_by,
                  ignore_errors)


    

    saveas_indices = {'ages': 0, 'fulltree': 1, 'subtrees': 2}
    saveas_i = saveas_indices[saveas]
    save_result = savefunctions[saveas]

    global showtree
    showtree = def_showtree(measures, show)

    phyltree = PhylTree.PhylogeneticTree(phyltreefile.format(ensembl_version))
    
    with Stream(outfile, 'w') as out:
        if saveas == 'ages':
            header = ['name'] + ['branch_'+m for m in measures] + \
                     ['age_'+m for m in measures] + \
                     ['calibrated', 'type', 'parent', 'taxon', 'root', 'subgenetree']
            out.write('\t'.join(header) + '\n')

        for i, mlcfile in enumerate(mlcfiles, start=1):
            percentage = float(i) / nb_mlc * 100
            print("\r%5d/%-5d (%3.2f%%) %s" % (i, nb_mlc, percentage, mlcfile),
                  end=' ')
            try:
                result = process(mlcfile, ensembl_version, phyltree,
                                 replace_nwk, replace_by, measures,
                                 method2=method2, unweighted=unweighted,
                                 tocalibrate=tocalibrate, keeproot=keeproot)
                save_result(result[saveas_i], out)
            except BaseException as err:
                print()
                if ignore_errors:
                    logger.error("Skip %r: %r", mlcfile, err)
                else:
                    raise
    print()
    

def readfromfiles(filenames):  # ~~> CLItools
    lines = []
    for filename in filenames:
        with open(filename) as f:
            lines += [line.rstrip() for line in f if not line.startswith('#')]
    return lines


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__, 
                                formatter_class=argparse.RawTextHelpFormatter)
    
    ### Input options
    gi = parser.add_argument_group('INPUT PARAMETERS')

    gi.add_argument('outfile')
    gi.add_argument('mlcfiles', nargs='+')
    gi.add_argument('--fromfile', action='store_true',
                    help='Take mlcfiles from a file (commented lines omitted).')
    gi.add_argument('-e', '--ensembl-version', type=int,
                    default=ENSEMBL_VERSION,
                    help='[%(default)s]')
    gi.add_argument('-p', '--phyltreefile', default=PHYLTREEFILE,
                    help='[%(default)s]')
    gi.add_argument('-r', '--replace-nwk', default='\.mlc$',
                    help='string to be replaced by REPLACE_BY to find the'\
                         ' tree file [%(default)s]')
    gi.add_argument('-R', '--replace-by', default='.nwk',
                    help='replacement string file [%(default)s]')
    
    ### Run options
    gr = parser.add_argument_group('RUN PARAMETERS')
    
    gr.add_argument('-2', '--method2', action='store_true')
    gr.add_argument('-u', '--unweighted', action='store_true', 
                    help='average weighted by the size of clusters')
    gr.add_argument('-c', '--tocalibrate', default="isdup", metavar='TESTS',
                    help=get_tocalibrate.__doc__ + '[%(default)r]')
    gr.add_argument('-k', '--keeproot', action='store_true', 
                    help="Calibrate data immediately following the root " \
                         "(not recommended).")
    gr.add_argument("-i", "--ignore-errors", action="store_true", 
                    help="On error, print the error and continue the loop.")
    
    ### Output options
    go = parser.add_argument_group('OUTPUT PARAMETERS')
    
    go.add_argument('--measures', nargs='*',
                    default=['t', 'dN', 'dS', 'dist'],
                    choices=['t', 'dN', 'dS', 'dist'],
                    help='Which distance measure: dist (from the newick ' \
                         'tree) or dS,dN,t (from codeml)')
    go.add_argument('-t', '--tofulltree', dest='saveas', action='store_const',
                    const='fulltree', default='ages',
                    help='Do not compute the table, but save trees in one'\
                        ' newick file with the chosen measure as distance')
    go.add_argument('-s', '--tosubtrees', dest='saveas', action='store_const',
                    const='subtrees', default='ages',
                    help='Do not compute the table, but save subtrees in '\
                        'one newick file with the chosen measure as ' \
                        'distance. "subtree" means the gene subtree '\
                        'contained between two speciations.')
    
    ### Display options
    gd = parser.add_argument_group('Display parameters')
    
    gd.add_argument('--show',
            help=('"gui": start the interactive ete3 tree browser\n'
                  '"notebook": create global var FULLTREE and TS for rendering\n'
                  'other value: save to file'))
    gd.add_argument('-v', '--verbose', action='store_true', 
                    help='print progression along tree')
    args = parser.parse_args()
    dictargs = vars(args)
    if dictargs.pop('fromfile'):
        dictargs['mlcfiles'] = readfromfiles(dictargs['mlcfiles'])
    main(**dictargs)

