#!/usr/bin/env python

"""Parse results of codeml (.mlc files) to save dN, dS values in a table"""

from __future__ import print_function

import sys
import os.path
import re
import argparse
import ete3
import LibsDyogen.utils.myPhylTree as PhylTree

from select_leaves_from_specieslist import convert_gene2species


def print_if_verbose(*args, **kwargs):
    pass

def showtree(fulltree, ages):
    pass

def def_showtree(show=False):
#    fulltree.show()
    if show:
        def showtree(fulltree, ages):
            ages_dict = {row[0]: row[1:] for row in ages}
            for n in fulltree.traverse():
                info = ages_dict.get(n.name, [0, 'leaf'])
                n.add_feature('age', info[0])
                n.add_feature('type', info[1])
            fulltree.show()
    else:
        def showtree(fulltree, ages):
            pass
    return showtree


def branch2nb(mlc, replace_nwk='.mlc'):
    """mlc: filehandle"""
    print_if_verbose()
    regex = re.compile(r'^(.*);$')
    regex_lnL = re.compile(r'^lnL\(')
    # get the line listing all branches
    line = mlc.next().rstrip()
    while not regex_lnL.match(line):
        line = mlc.next().rstrip()
    branches_line = mlc.next().rstrip()
    branches = [b.split('..') for b in branches_line.split()]
    # get translation by looking at the newick tree lines.
    # I could also get it by the lines #1: ENS...
    line = mlc.next().rstrip()
    while not regex.match(line):
        line = mlc.next().rstrip()
    tree_nbs = ete3.Tree(line)
    line = mlc.next().rstrip()
    while not regex.match(line):
        line = mlc.next().rstrip()
    tree_ids = ete3.Tree(line)

    id2nb = dict(zip(tree_ids.get_leaf_names(), tree_nbs.get_leaf_names()))
    nb2id = dict(zip(tree_nbs.get_leaf_names(), tree_ids.get_leaf_names()))

    # get internal nodes nb-to-id conversion
    nwkfile = mlc.name.replace(replace_nwk, '.nwk')
    try:
        fulltree = ete3.Tree(nwkfile, format=1)
    except ete3.parser.newick.NewickError as e:
        if os.path.exists(nwkfile):
            print("\nNewickError: Malformed newick tree structure in %r" \
                    % nwkfile,
                    file=sys.stderr)
        else:
            print("\nNewickError: Unexisting tree file %r" % nwkfile,
                    file=sys.stderr)
        sys.exit(1)


    for leafnb, leafid in zip(tree_nbs.get_leaves(), fulltree.get_leaves()):
        leafid.add_feature('nb', leafnb.name)

    # branches follow the order of the newick string.
    while branches:
        base, tip = branches.pop()
        print_if_verbose("%-8s" % (base + '..' + tip), end=' ')
        try:
            tip_id = nb2id[tip]
            tip_id_anc = fulltree.search_nodes(name=tip_id)[0].get_ancestors()
            try:
                base_id = tip_id_anc[0].name
                nb2id[base] = base_id
                # Add number in the fulltree:
                tip_id_anc[0].add_feature('nb', base)
                # Also update the tree_nbs
                base_nb_node = tree_nbs.search_nodes(name=tip)[0].get_ancestors()[0]
                base_nb_node.name = base
                print_if_verbose('ok')
            except IndexError:
                print('root (%r: %r) cannot have ancestors' % (tip, tip_id))
                pass
        except KeyError as e:
            branches = [base, tip] + branches
            print_if_verbose('KeyError')

    print_if_verbose(tree_nbs.get_ascii())
    return nb2id, tree_nbs, fulltree


def get_dNdS(mlc):
    """mlc: filehandle"""
    #reg_dNdS = re.compile(r'dN & dS for each branch')
    reg_dNdS = re.compile(r'\s+branch\s+t\s+N\s+S\s+dN/dS\s+dN\s+dS\s+N\*dN\s+S\*dS')
    line = mlc.next().rstrip()
    while not reg_dNdS.match(line):
        line = mlc.next().rstrip()
    assert mlc.next().rstrip() == ''  # skip blank line
    dNdS = {}
    line = mlc.next().rstrip()
    while line != '':
        #print line
        fields = line.split()
        dNdS[fields[0]] = [float(x) for x in fields[1:]]
        line = mlc.next().rstrip()
    return dNdS
    

def sum_average_dNdS(dNdS, nb2id, tree_nbs):
    """Compute the average dS depth of each node. (between this node and leaves)"""
    dS, dN = {}, {}
    for n in tree_nbs.traverse('postorder'):
        node = n.name
        if n.is_leaf():
            dS[nb2id[node]] = 0
            dN[nb2id[node]] = 0
            dS[node] = 0
            dN[node] = 0
        else:
            dS_children = [dS[nb2id[c.name]] + dNdS[node + '..' + c.name][5] for c in n.children]
            dN_children = [dN[nb2id[c.name]] + dNdS[node + '..' + c.name][4] for c in n.children]
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
    ancgene2sp = re.compile(r'([A-Z][A-Za-z_.-]+)ENSGT')
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
                taxon = ancgene2sp.match(scname).group(1).replace('.', ' ')
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
                ages.append([scname, str(node_age), "spe"])
                # compute average dS to each next speciation.
                # update age of each posterior dup
                for c in n.children:
                    # walk descendants until speciation.
                    # need to get the age of next speciation and compute the time
                    # between the two speciation.
                    # NOTE: could actually have used the phyltree.branches
                    branch_length = node_age - subtree[c.name]['age']
                    scaling_dS = subtree[c.name]['dS'] + \
                                    dNdS[node + '..' + c.name][5]
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
                            ages.append([nb2id[nextnode.name], str(age), "dup"])
                            nextnodes.extend(nextnode.children)
                        subtree.pop(nextnode.name)

                # then reset dS to zero
                subtree[node]['dS'] = 0
    return ages


def bound_average_dS_2(dNdS, nb2id, tree_nbs, phyltree):
    """each dS of a duplication is normalized by the D' + dS, where:
        - dS: average of branch length leading to next speciation, in dS
        - D': branch length from the previous speciation to this duplication"""
    pass


def save_ages(ages, opened_outfile):
    for row in ages:
        opened_outfile.write("\t".join(row) + "\n")


def process_mlc(mlcfile, phyltree, replace_nwk='.mlc'):
    with open(mlcfile) as mlc:
        nb2id, tree_nbs, fulltree = branch2nb(mlc, replace_nwk)
        dNdS = get_dNdS(mlc)

    #dN, dS = sum_average_dNdS(dNdS, nb2id, tree_nbs)
    ages = bound_average_dS(dNdS, nb2id, tree_nbs, phyltree)
    showtree(fulltree, ages)
    return ages


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('outfile')
    parser.add_argument('mlcfiles', nargs='+')
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='print progression along tree')
    parser.add_argument('--show', action='store_true',
                        help='start the interactive ete3 tree browser')
    parser.add_argument('-r', '--replace-nwk', default='.mlc',
                        help='string to be replaced by .nwk to find the tree'\
                               ' file [%(default)s]')
    args = parser.parse_args()
    outfile = args.outfile
    mlcfiles = args.mlcfiles
    nb_mlc = len(mlcfiles)
    
    # "compile" some functions to avoid redondant tests ("if verbose: ...")
    if args.verbose:
        def print_if_verbose(*args, **kwargs):
            print(*args, **kwargs)

    showtree = def_showtree(args.show)

    phyltree = PhylTree.PhylogeneticTree("/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data85/PhylTree.Ensembl.85.conf")
    
    with open(outfile, 'a') as out:
        for i, mlcfile in enumerate(mlcfiles, start=1):
            percentage = float(i) / nb_mlc * 100
            print("\r%5d/%-5d (%3.2f%%) %s" % (i, nb_mlc, percentage, mlcfile),
                    file=sys.stderr, end=' ')
            try:
                ages = process_mlc(mlcfile, phyltree, args.replace_nwk)
                save_ages(ages, out)
            except:
                print()
                raise
    print()
    
