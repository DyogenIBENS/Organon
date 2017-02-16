#!/usr/bin/env python3

"""Parse results of codeml (.mlc files) to save dN, dS values in a table"""

from __future__ import print_function

import sys
import os.path
import re
import argparse
import ete3
import LibsDyogen.myPhylTree as PhylTree # my custom python3 version

from select_leaves_from_specieslist import convert_gene2species


ANCGENE2SP = re.compile(r'([A-Z][A-Za-z_.-]+)ENS')


def print_if_verbose(*args, **kwargs):
    """Default print function. Assume that --verbose is False, so print nothing."""
    pass

def showtree(fulltree, ages):
    """Default `showtree` function: do nothing. This function can be
    redefined using `def_showtree`"""
    pass

def def_showtree(show=False):
    """Depending on the boolean argument 'show', return the function `showtree`:
        if show=False, this function does nothing;
        else, showtree display an interactive tree using ete3."""
    if show:
        def showtree(fulltree, ages):
            ages_dict = {row[0]: row[1:] for row in ages}

            # define duplication node style:
            ns_dup = ete3.NodeStyle()
            ns_dup['fgcolor'] = 'red'
            
            for n in fulltree.traverse():
                info = ages_dict.get(n.name, [0.0, 'leaf'])
                n.add_feature('age', info[0])
                n.add_feature('type', info[1])
                if info[1] == 'dup':
                    n.set_style(ns_dup)
                # Add elements to be rendered
                #if not n.is_leaf():
                #    n.add_face(ete3.TextFace(n.name), column=0)
                n.add_face(ete3.AttrFace('age',
                                         text_prefix='age: ',
                                         fgcolor='grey',
                                         formatter='%.4g'),
                           column=0,
                           position='branch-bottom')
                if hasattr(n, 'dS'):
                    n.add_face(ete3.AttrFace('dS',
                                             text_prefix='dS = ',
                                             fgcolor='grey',
                                             formatter='%.4g'),
                               column=0,
                               position='branch-bottom')
            # define a tree style:
            ts = ete3.TreeStyle()
            ts.show_branch_length = True

            fulltree.show(tree_style=ts)
    else:
        def showtree(fulltree, ages):
            pass
    return showtree


def load_fulltree(mlcfile, replace_nwk='.mlc'):
    """return the ete3.Tree object corresponding to a given .mlc file.
    Catch errors (file does not exist / wrong format)"""
    nwkfile = mlcfile.replace(replace_nwk, '.nwk')
    try:
        return ete3.Tree(nwkfile, format=1)
    except ete3.parser.newick.NewickError as e:
        if os.path.exists(nwkfile):
            print("\nNewickError: Malformed newick tree structure in %r" \
                    % nwkfile,
                    file=sys.stderr)
        else:
            print("\nNewickError: Unexisting tree file %r" % nwkfile,
                    file=sys.stderr)
        sys.exit(1)


def branch2nb(mlc, replace_nwk='.mlc'):
    """Parse the codeml result file (.mlc) to return 2 trees:
    tree_nbs: the tree with node labelled as numbers.
    tree_ids: the tree with original node labels (including inner nodes).
    
    Also return the dictionaries to convert nb -> id and conversely.

    Arguments:
        - mlc: an opened file"""
    print_if_verbose()
    regex = re.compile(r'^(.*);$')
    regex_lnL = re.compile(r'^lnL\(')
    # get the line listing all branches
    line = mlc.readline().rstrip()
    while not regex_lnL.match(line):
        line = mlc.readline().rstrip()
    branches_line = mlc.readline().rstrip()
    branches = [b.split('..') for b in branches_line.split()]
    # get translation by looking at the newick tree lines.
    # I could also get it by the lines #1: ENS...
    line = mlc.readline().rstrip()
    while not regex.match(line):
        line = mlc.readline().rstrip()
    tree_nbs = ete3.Tree(line)

    line = mlc.readline().rstrip()
    while not regex.match(line):
        line = mlc.readline().rstrip()
    tree_ids = ete3.Tree(line)

    id2nb = dict(zip(tree_ids.get_leaf_names(), tree_nbs.get_leaf_names()))
    nb2id = dict(zip(tree_nbs.get_leaf_names(), tree_ids.get_leaf_names()))

    # get internal nodes nb-to-id conversion (need the tree with internal node
    # annotations)
    fulltree = load_fulltree(mlc.name, replace_nwk)

    for leafnb, leafid in zip(tree_nbs.get_leaves(), fulltree.get_leaves()):
        leafid.add_feature('nb', leafnb.name)

    # branches follow the order of the newick string.
    while branches:
        base, tip = branches.pop()
        print_if_verbose("%-8s" % (base + '..' + tip), end=' ')
        try:
            tip_id = nb2id[tip]
            try:
                base_node = fulltree.search_nodes(name=tip_id)[0].up
                while len(base_node.children) == 1:
                    base_node = base_node.up
                base_id = base_node.name
                print_if_verbose("%s -> %s  " % (base_id, tip_id), end=" ")
                nb2id[base] = base_id
                id2nb[base_id] = base
                # Add number in the fulltree:
                base_node.add_feature('nb', base)
                # Also update the tree_nbs
                base_nb_node = tree_nbs.search_nodes(name=tip)[0].up
                base_nb_node.name = base
                print_if_verbose('ok')
            except AttributeError:
                print('root (%r: %r) cannot have ancestors' % (tip, tip_id))
                pass
        except KeyError as e:
            # Not found now, put it back in the queue for later
            branches.insert(0, (base, tip))
            print_if_verbose('KeyError')

    print_if_verbose(tree_nbs.get_ascii())
    return id2nb, nb2id, tree_nbs, fulltree


def get_dNdS(mlc):
    """Parse table of dN/dS from codeml output file.
    
    mlc: filehandle
    """
    #reg_dNdS = re.compile(r'dN & dS for each branch')
    reg_dNdS = re.compile(r'\s+branch\s+t\s+N\s+S\s+dN/dS\s+dN\s+dS\s+N\*dN\s+S\*dS')
    line = mlc.readline().rstrip()
    while not reg_dNdS.match(line):
        line = mlc.readline().rstrip()
    assert mlc.readline().rstrip() == ''  # skip blank line
    dNdS = {}
    line = mlc.readline().rstrip()
    while line != '':
        #print line
        fields = line.split()
        dNdS[fields[0]] = [float(x) for x in fields[1:]]
        line = mlc.readline().rstrip()
    return dNdS
    

def set_dS_fulltree(fulltree, id2nb, dNdS):
    """Add codeml dS on each branch of the complete tree (with missing species
    nodes).
    When a node has a single child: take the branch dS on which it is, 
    and divide proportionnally to each segment dist."""
    for node in fulltree.get_descendants('postorder'): # exclude root
        if len(node.children) != 1:
            parent = node.up
            intermediates = []
            while parent.up and len(parent.children) == 1: # under root, up is None
                intermediates.append(parent)
                parent = parent.up

            if parent.is_root():
                # then this is a root with a single children: cannot have dS
                dS_tot = 0
            else:
                dS_tot = dNdS[id2nb[parent.name] + '..' + id2nb[node.name]][5]
    
            if intermediates: # this if is not necessary, but maybe for efficiency
                dist_tot = sum(n.dist for n in intermediates + [node])

                if dist_tot == 0:
                    print("WARNING: dS = 0 between %r and %r" % (parent.name,
                            node.name),
                            file=sys.stderr)
                    dist_tot = 1. # should set to NaN to be sure.

                for inter_node in intermediates + [node]:
                    inter_node.add_feature('dS', dS_tot * inter_node.dist / dist_tot)
            else:
                node.add_feature('dS', dS_tot)


def rm_erroneous_ancestors(fulltree, phyltree):
    """Some ancestors with only one child are present in the species phylogeny.
    They must be remomved when they have an age of zero"""
    for node in fulltree.iter_descendants():
        if not node.is_leaf():
            try:
                taxon = ANCGENE2SP.match(node.name).group(1).replace('.', ' ')
            except AttributeError:
                raise RuntimeError("Can not match species name in %r" % node.name)
            if phyltree.ages[taxon] == 0 and hasattr(node, 'reinserted'):
                node.delete(prevent_nondicotomic=False, preserve_branch_length=True)



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


def bound_average_dS(dNdS, nb2id, tree_nbs, phyltree):
    """calibrate tree using speciation nodes."""
    ages = []
    subtree = {}
    for n in tree_nbs.traverse('postorder'):
        node = n.name
        scname = nb2id[node] # scientific name
        print_if_verbose("%2s. %s:" % (node, scname), end=' ')
        if n.is_leaf():
            taxon = convert_gene2species(scname)
            subtree[node] = {'taxon': taxon, 'dS': 0, 'age': 0} #, 'next_age': 0}
            #ages[scname] = 0
            print_if_verbose("Leaf")
        else:
            try:
                taxon = ANCGENE2SP.match(scname).group(1).replace('.', ' ')
            except AttributeError:
                raise RuntimeError("Can not match species name in %r" % scname)

            subtree[node] = {'taxon': taxon}
            # determine if speciation or not
            children_taxa = set(subtree[c.name]['taxon'] for c in n.children)
            # compute and store average dS to next spe, for each child
            children_dS = [subtree[c.name]['dS'] + \
                            dNdS[node + '..' + c.name][5] \
                            for c in n.children]
            if len(children_taxa) == 1:
                # It is a duplication: average children_dS
                subtree[node]['dS'] = sum(children_dS) / len(children_dS)
                print_if_verbose("Dupl; dS=%s" % subtree[node]['dS'])
                # store the age of the next speciation event.
                # Since it is a duplication, this should be the same for
                # children[0] and children[1].
                # This 'age' will later be modified (rescaled according to dS)
                subtree[node]['age'] = subtree[n.children[0].name]['age']
            else:
                print_if_verbose("Spe")
                # speciation: store the age of this taxon
                node_age = subtree[node]['age'] = phyltree.ages[taxon]
                ages.append([scname, node_age, "spe"])
                # compute average dS to each next speciation.
                # update age of each posterior dup
                if len(n.children) > 1:
                    pass
                for c in n.children:
                    # walk descendants until speciation.
                    # need to get the age of next speciation and compute the time
                    # between the two speciation.
                    # NOTE: could actually have used the phyltree.branches
                    branch_length = node_age - subtree[c.name]['age']
                    dS_genebranch = dNdS[node + '..' + c.name][5]
                    scaling_dS = subtree[c.name]['dS'] + dS_genebranch
                    print_if_verbose("    climb up to next speciation: " \
                                      "scaling_dS=%s; br_len=%s" % \
                                        (scaling_dS, branch_length))
                    nextnodes = [c]
                    while nextnodes:
                        nextnode = nextnodes.pop(0)
                        nextnode_dS = subtree[nextnode.name]['dS']
                        print_if_verbose("    - %2s. %s: dS=%s" % \
                                            (nextnode.name,
                                             nb2id[nextnode.name],
                                             nextnode_dS))
                        if nextnode_dS > 0:
                            age = node_age - (1 - nextnode_dS / scaling_dS) * branch_length
                            #ages[nextnode.name] = age
                            ages.append([nb2id[nextnode.name], age, "dup"])
                            nextnodes.extend(nextnode.children)
                        subtree.pop(nextnode.name)

                # then reset dS to zero
                subtree[node]['dS'] = 0
    return ages


def bound_average_dS(dNdS, id2nb, fulltree, phyltree):
    ages = []
    subtree = {}
    #for node in fulltree.traverse('postorder'):
    #    print(id2nb.get(node.name), end=' ')
    #print()
    #fulltree.show()
    for node in fulltree.traverse('postorder'):
        scname = node.name
        print_if_verbose("* %3s. %s:" % (id2nb.get(scname), scname), end=' ')
        if node.is_leaf():
            taxon = convert_gene2species(scname)
            subtree[scname] = {'taxon': taxon, 'tmp_dS': 0, 'age': 0} #, 'next_age': 0}
            #ages[scname] = 0
            print_if_verbose("Leaf")
        else:
            try:
                taxon = ANCGENE2SP.match(scname).group(1).replace('.', ' ')
            except AttributeError:
                raise RuntimeError("Can not match species name in %r" % scname)

            subtree[scname] = {'taxon': taxon}
            children_taxa = set((subtree[ch.name]['taxon'] for ch in node.children))

            if len(children_taxa & set((taxon,))) == 1:
                # it is a duplication:
                children_dS = [ch.dS + subtree[ch.name]['tmp_dS'] for ch in node.children]
                subtree[scname]['tmp_dS'] = sum(children_dS) / len(children_dS)

                print_if_verbose("Dupl; dS=%s" % subtree[scname]['tmp_dS'])
                # store the age of the next speciation event.
                # Since it is a duplication, this should be the same for
                # children[0] and children[1].
                # This 'age' will later be modified (rescaled according to dS)
                subtree[scname]['age'] = subtree[node.children[0].name]['age']

            else:
                # it is a speciation.
                print_if_verbose("Spe")
                # store the age of this taxon
                node_age = subtree[scname]['age'] = phyltree.ages[taxon]
                ages.append([scname, node_age, "spe"])
                # climb up tree and assign an age to each duplication
                for ch in node.children:
                    # walk descendants until speciation.
                    # need to get the age of next speciation and compute the time
                    # between the two speciation.
                    # NOTE: could actually have used the phyltree.branches
                    branch_length = node_age - subtree[ch.name]['age']
                    dS_genebranch = ch.dS
                    scaling_dS = subtree[ch.name]['tmp_dS'] + dS_genebranch
                    print_if_verbose("    climb up to next speciation: " \
                                      "scaling_dS=%s; br_len=%s" % \
                                        (scaling_dS, branch_length))
                    nextnodes = [ch]
                    while nextnodes:
                        nextnode = nextnodes.pop(0)
                        try:
                            nextnode_dS = subtree[nextnode.name]['tmp_dS']
                        except KeyError as err:
                            err.args += ("Error: Node exists twice in the tree.",
                                         "You may need to rerun `prune2family.py`",)
                            raise
                        print_if_verbose("    - %2s. %s: dS=%s" % \
                                            (id2nb.get(nextnode.name),
                                             nextnode.name,
                                             nextnode_dS))
                        if nextnode_dS > 0:
                            age = node_age - (1 - nextnode_dS / scaling_dS) * branch_length
                            ages.append([nextnode.name, age, "dup"])
                            nextnodes.extend(nextnode.children)
                        #print("    Pop: %s" % nextnode.name)
                        #print("    nextnodes: ", nextnodes)
                        #print("    next children: ", nextnode.children)
                        subtree.pop(nextnode.name)

                # then reset dS to zero
                subtree[scname]['tmp_dS'] = 0
    return ages


def bound_average_dS_2(dNdS, id2nb, fulltree, phyltree):
    """each dS of a duplication is normalized by the D' + dS, where:
        - dS: average of branch length leading to next speciation, in dS
        - D': branch length from the previous speciation to this duplication"""
    ages = []
    subtree = {}
    for node in fulltree.traverse('postorder'):
        scname = node.name
        print_if_verbose("* %3s. %s:" % (id2nb.get(scname), scname), end=' ')
        if node.is_leaf():
            taxon = convert_gene2species(scname)
            subtree[scname] = {'taxon': taxon, 'tmp_dS': 0, 'age': 0} #, 'next_age': 0}
            #ages[scname] = 0
            print_if_verbose("Leaf")
        else:
            try:
                taxon = ANCGENE2SP.match(scname).group(1).replace('.', ' ')
            except AttributeError:
                raise RuntimeError("Can not match species name in %r" % scname)

            subtree[scname] = {'taxon': taxon}
            #for child in node.children:
            #    child_name = child.name
            #    child_taxon = subtree[child_name]
            children_taxa = set((subtree[ch.name]['taxon'] for ch in node.children))

            if len(children_taxa & set((taxon,))) == 1:
                # it is a duplication:
                children_dS = [ch.dS + subtree[ch.name]['tmp_dS'] for ch in node.children]
                subtree[scname]['tmp_dS'] = sum(children_dS) / len(children_dS)

                print_if_verbose("Dupl; dS=%s" % subtree[scname]['tmp_dS'])
                # store the age of the next speciation event.
                # Since it is a duplication, this should be the same for
                # children[0] and children[1].
                # This 'age' will later be modified (rescaled according to dS)
                subtree[scname]['age'] = subtree[node.children[0].name]['age']

            else:
                # it is a speciation.
                print_if_verbose("Spe")
                # store the age of this taxon
                node_age = subtree[scname]['age'] = phyltree.ages[taxon]
                ages.append([scname, node_age, "spe"])
                # climb up tree and assign an age to each duplication
                for ch in node.children:
                    ### This is where version _2 is different
                    # walk descendants until speciation.
                    # need to get the age of next speciation and compute the time
                    # between the two speciation.
                    branch_length = node_age - subtree[ch.name]['age']
                    print_if_verbose("    climb up to next speciation: " \
                                      "br_len=%s" % branch_length)
                    nextnodes = [(ch, ch.dS)]
                    while nextnodes:
                        nextnode, next_path_dS = nextnodes.pop(0)
                        nextnode_dS = subtree[nextnode.name]['tmp_dS']
                        print_if_verbose("    - %2s. %s: dS(to speciation)=%s"\
                                         "; dS(from speciation)=%s" % \
                                            (id2nb.get(nextnode.name),
                                             nextnode.name,
                                             nextnode_dS, next_path_dS))
                        if nextnode_dS > 0:
                            scaling_dS = next_path_dS + nextnode_dS
                            age = node_age - (1 - nextnode_dS / scaling_dS) * branch_length
                            #ages[nextnode.name] = age
                            ages.append([nextnode.name, age, "dup"])
                            nextnodes.extend((nch, nch.dS + next_path_dS) \
                                                for nch in nextnode.children)
                        subtree.pop(nextnode.name)

                # then reset dS to zero
                subtree[scname]['tmp_dS'] = 0
    return ages


def bound_average_2(fulltree, phyltree, measure='dS'):
    """normalize duplication node position between speciation nodes.
     scaling: d / (D' + d)
        - d: average of branch length leading to next speciation, in 'measure'.
             'measure' must be an attribute present in fulltree ('dS', 'dist').
        - D': branch length from the previous speciation to this duplication"""
    ages = []
    subtree = {} # temporary subtree while traversing from one speciation to another
    for node in fulltree.traverse('postorder'):
        scname = node.name
        print_if_verbose("* %s:" % scname, end=' ')
        if node.is_leaf():
            taxon = convert_gene2species(scname)
            subtree[scname] = {'taxon': taxon, 'tmp_m': 0, 'age': 0} #, 'next_age': 0}
            #ages[scname] = 0
            print_if_verbose("Leaf")
        else:
            try:
                taxon = ANCGENE2SP.match(scname).group(1).replace('.', ' ')
            except AttributeError:
                raise RuntimeError("Can not match species name in %r" % scname)

            subtree[scname] = {'taxon': taxon}
            children_taxa = set((subtree[ch.name]['taxon'] for ch in node.children))

            if len(children_taxa & set((taxon,))) == 1:
                # it is a duplication:
                children_m = [ch.__getattribute__(measure) +
                              subtree[ch.name]['tmp_m']
                                  for ch in node.children]
                subtree[scname]['tmp_m'] = sum(children_m) / len(children_m)

                print_if_verbose("Dupl; m=%s" % subtree[scname]['tmp_m'])
                # store the age of the next speciation event.
                # Since it is a duplication, this should be the same for
                # children[0] and children[1].
                # This 'age' will later be modified (rescaled according to dS)
                subtree[scname]['age'] = subtree[node.children[0].name]['age']

            else:
                # it is a speciation.
                print_if_verbose("Spe")
                # store the age of this taxon
                node_age = subtree[scname]['age'] = phyltree.ages[taxon]
                ages.append([scname, node_age, "spe"])
                # climb up tree and assign an age to each duplication
                for ch in node.children:
                    ### This is where version _2 is different
                    # walk descendants until speciation.
                    # need to get the age of next speciation and compute the time
                    # between the two speciation.
                    branch_length = node_age - subtree[ch.name]['age']
                    print_if_verbose("    climb up to next speciation: " \
                                      "br_len=%s" % branch_length)
                    nextnodes = [(ch, ch.__getattribute__(measure))]
                    while nextnodes:
                        nextnode, next_path_m = nextnodes.pop(0)
                        try:
                            nextnode_m = subtree[nextnode.name]['tmp_m']
                        except KeyError as err:
                            err.args += ("Error: Node exists twice in the tree.",
                                         "You may need to rerun `prune2family.py`",)
                            raise
                        print_if_verbose("    - %s: measure(to speciation)=%s"\
                                         "; measure(from speciation)=%s" % \
                                            (nextnode.name,
                                             nextnode_m, next_path_m))
                        if nextnode_m > 0:
                            scaling_m = next_path_m + nextnode_m
                            age = node_age - (1 - nextnode_m / scaling_m) * branch_length
                            #ages[nextnode.name] = age
                            ages.append([nextnode.name, age, "dup"])
                            nextnodes.extend((nch,
                                              nch.__getattribute__(measure) + next_path_m) \
                                                for nch in nextnode.children)
                        subtree.pop(nextnode.name)

                # then reset dS to zero
                subtree[scname]['tmp_m'] = 0
    return ages


def save_ages(ages, opened_outfile):
    for row in ages:
        row[1] = str(row[1])
        opened_outfile.write("\t".join(row) + "\n")


def process_mlc(mlcfile, phyltree, replace_nwk='.mlc'):
    """main command: scale ages based on *dS* (method1)"""
    with open(mlcfile) as mlc:
        id2nb, nb2id, tree_nbs, fulltree = branch2nb(mlc, replace_nwk)
        dNdS = get_dNdS(mlc)

    #dN, dS = sum_average_dNdS(dNdS, nb2id, tree_nbs)
    rm_erroneous_ancestors(fulltree, phyltree)
    set_dS_fulltree(fulltree, id2nb, dNdS)
    ages = bound_average_dS(dNdS, id2nb, fulltree, phyltree)
    showtree(fulltree, ages)
    return ages

def process_mlc2(mlcfile, phyltree, replace_nwk='.mlc'):
    """main command: scale ages based on *dS* (method2)"""
    with open(mlcfile) as mlc:
        id2nb, nb2id, tree_nbs, fulltree = branch2nb(mlc, replace_nwk)
        dNdS = get_dNdS(mlc)

    #dN, dS = sum_average_dNdS(dNdS, nb2id, tree_nbs)
    rm_erroneous_ancestors(fulltree, phyltree)
    set_dS_fulltree(fulltree, id2nb, dNdS)
    #ages = bound_average_dS_2(dNdS, id2nb, fulltree, phyltree)
    ages = bound_average_2(fulltree, phyltree, measure='dS')
    showtree(fulltree, ages)
    return ages


def process_nwk(mlcfile, phyltree, replace_nwk='.mlc'):
    """main command: scale ages based on *dist* (method1)"""
    raise NotImplementedError


def process_nwk2(mlcfile, phyltree, replace_nwk='.mlc'):
    """main command: scale ages based on *dist* (method1)"""
    fulltree = load_fulltree(mlcfile, replace_nwk)
    rm_erroneous_ancestors(fulltree, phyltree)
    ages = bound_average_2(fulltree, phyltree, measure='dist')
    showtree(fulltree, ages)
    return ages


class Out(object):
    """Context Manager class (for use with `with` statement). Do the exact
    same as `open()`, but if filename is '-', open stdout for writing."""
    write_modes = ('w', 'x', 'a')

    def __init__(self, filename, *args, **kwargs):
        self.filename = filename
        if filename == '-':
            mode = args[0] if args else kwargs.get('mode', 'r')
            if any(letter in mode for letter in write_modes):
                self.file_obj = sys.stdout
            else:
                self.file_obj = sys.stdin
        else:
            self.file_obj = open(filename, *args, **kwargs)

    def __enter__(self):
        return self.file_obj

    def __exit__(self, type, value, traceback):
        if self.filename != '-':
            self.file_obj.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('outfile')
    parser.add_argument('mlcfiles', nargs='+')
    parser.add_argument('--method2', action='store_true')
    parser.add_argument('--measure', default='dS', choices=['dS', 'dist'],
                        help='which distance measure: dS (from codeml) or ' \
                             'dist (from PhyML)')
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='print progression along tree')
    parser.add_argument('--show', action='store_true',
                        help='start the interactive ete3 tree browser')
    parser.add_argument('-r', '--replace-nwk', default='.mlc',
                        help='string to be replaced by .nwk to find the tree'\
                               ' file [%(default)s]')
    parser.add_argument("-i", "--ignore-errors", action="store_true", 
                        help="On error, print the error and continue the loop.")
    args = parser.parse_args()
    outfile = args.outfile
    mlcfiles = args.mlcfiles
    nb_mlc = len(mlcfiles)
    
    # "compile" some functions to avoid redondant tests ("if verbose: ...")
    if args.verbose:
        def print_if_verbose(*args, **kwargs):
            print(*args, **kwargs)

    if args.measure == 'dS':
        process = process_mlc2 if args.method2 else process_mlc
    elif args.measure == 'dist':
        process = process_nwk2 if args.method2 else process_nwk
    else:
        raise RuntimeError("Unknown measure %r" % args.measure)

    showtree = def_showtree(args.show)

    phyltree = PhylTree.PhylogeneticTree("/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data85/PhylTree.Ensembl.85.conf")
    
    with Out(outfile, 'w') as out:
        for i, mlcfile in enumerate(mlcfiles, start=1):
            percentage = float(i) / nb_mlc * 100
            print("\r%5d/%-5d (%3.2f%%) %s" % (i, nb_mlc, percentage, mlcfile),
                  end=' ')
            try:
                ages = process(mlcfile, phyltree, args.replace_nwk)
                save_ages(ages, out)
            except BaseException as err:
                print()
                if args.ignore_errors:
                    print("Skip %r: %r" % (mlcfile, err), file=sys.stderr)
                else:
                    raise
    print()
    
