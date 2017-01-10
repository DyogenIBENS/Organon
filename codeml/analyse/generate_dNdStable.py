#!/usr/bin/env python2.7

"""Parse results of codeml (.mlc files) to save dN, dS values in a table"""

from __future__ import print_function

import sys
import re
import ete3
import LibsDyogen.utils.myPhylTree as PhylTree

from select_leaves_from_specieslist import convert_gene2species


def branch2nb(mlc):
    """mlc: filehandle"""
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
    nwkfile = mlc.name.replace('.mlc', '.nwk')
    fulltree = ete3.Tree(nwkfile, format=1)

    # branches follow the order of the newick string.
    while branches:
        base, tip = branches.pop()
        #print base + '..' + tip,
        try:
            tip_id = nb2id[tip]
            tip_id_anc = fulltree.search_nodes(name=tip_id)[0].get_ancestors()
            try:
                base_id = tip_id_anc[0].name
                nb2id[base] = base_id
                # Also update the tree_nbs
                base_nb_node = tree_nbs.search_nodes(name=tip)[0].get_ancestors()[0]
                base_nb_node.name = base
                #print 'ok'
            except IndexError:
                print('root (%r: %r) cannot have ancestors' % (tip, tip_id))
                pass
        except KeyError as e:
            branches = [base, tip] + branches
            #print 'KeyError'

    return nb2id, tree_nbs


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


def bound_average_dNdS(dNdS, nb2id, tree_nbs, phyltree):
    """calibrate tree using speciation nodes."""
    ancgene2sp = re.compile(r'([A-Z][a-z.]+)ENSGT')
    ages = {}
    subtree = {}
    for n in tree_nbs.traverse('postorder'):
        node = n.name
        scname = nb2id[node] # scientific name
        #print "%2s. %s:" % (node, scname),
        if n.is_leaf():
            taxon = convert_gene2species(scname)
            subtree[node] = {'taxon': taxon, 'dS': 0, 'age': 0} #, 'next_age': 0}
            ages[scname] = 0
            #print "Leaf"
        else:
            taxon = ancgene2sp.match(scname).group(1).replace('.', ' ')
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
                #print "Dupl; dS=%s" % subtree[node]['dS']
                # store the age of the next speciation event.
                # Since it is a duplication, this should be the same for
                # children[0] and children[1].
                # This 'age' will later be modified (rescaled according to dS)
                subtree[node]['age'] = subtree[n.children[0].name]['age']
            else:
                #print "Spe"
                # speciation: store the age of this taxon
                node_age = phyltree.ages[taxon]
                ages[scname] = subtree[node]['age'] = node_age
                # compute average dS to each next speciation.
                # update age of each posterior dup
                for c in n.children:
                    # walk descendants until speciation.
                    # need to get the age of next speciation and compute the time
                    # between the two speciation.
                    # NOTE: could actually have used the phyltree.branches
                    branch_length = node_age - subtree[c.name]['age']
                    scaling_dS = subtree[c.name]['dS'] + dNdS[node + '..' + c.name][5]
                    #print "    climb up to next speciation: scaling_dS=%s; br_len=%s" %(scaling_dS, branch_length)
                    nextnodes = [c]
                    while nextnodes:
                        nextnode = nextnodes.pop(0)
                        nextnode_dS = subtree[nextnode.name]['dS']
                        #print "    - %2s. %s: dS=%s" % (nextnode.name, nb2id[nextnode.name], nextnode_dS)
                        if nextnode_dS > 0:
                            age = node_age - (1 - nextnode_dS / scaling_dS) * branch_length
                            #ages[nextnode.name] = age
                            ages[nb2id[nextnode.name]] = age
                            nextnodes.extend(nextnode.children)
                        subtree.pop(nextnode.name)

                # then reset dS to zero
                subtree[node]['dS'] = 0
    return ages


def save_ages(ages, opened_outfile):
    for node, age in ages.iteritems():
        opened_outfile.write("%s\t%s\n" % (node, age))


def process_mlc(mlcfile, phyltree):
    with open(mlcfile) as mlc:
        nb2id, tree_nbs = branch2nb(mlc)
        dNdS = get_dNdS(mlc)

    #dN, dS = sum_average_dNdS(dNdS, nb2id, tree_nbs)
    ages = bound_average_dNdS(dNdS, nb2id, tree_nbs, phyltree)
    return ages


if __name__=='__main__':
    if len(sys.argv) < 3:
        print("USAGE: ./generate_dNdStable.py <outfile> <mlcfile> [<mlcfile>, ...]")
        sys.exit(1)

    outfile = sys.argv[1]
    phyltree = PhylTree.PhylogeneticTree("/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data85/PhylTree.Ensembl.85.conf")
    mlcfiles = sys.argv[2:]
    nb_mlc = len(mlcfiles)
    with open(outfile, 'a') as out:
        for i, mlcfile in enumerate(mlcfiles, start=1):
            percentage = float(i) / nb_mlc * 100
            print("%5d/%-5d (%3.2f%%) %s\r" % (i, nb_mlc, percentage, mlcfile),
                    file=sys.stderr, end='')
            ages = process_mlc(mlcfile, phyltree)
            save_ages(ages, out)
    print()
    
