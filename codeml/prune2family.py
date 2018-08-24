#!/usr/bin/env python3

from __future__ import print_function

"""Extract all gene subtrees rooted at a given ancestral taxon.
Also add missing speciation nodes."""


import re
import sys
import os.path
try:
    import argparse_custom as argparse
    #print('ARGPARSE CUSTOM')
except ImportError:
    import argparse
import multiprocessing as mp

from copy import copy

import ete3
import LibsDyogen.myPhylTree as PhylTree

# The 3 following imports are just so messy. TODO: write a unique conversion
# function, and/or centralize these functions in a single script.
from codeml.select_leaves_from_specieslist import convert_gene2species
from prot2gene import convert_prot2species
from seqtools.specify import load_conversion


ENSEMBL_VERSION = 85
PHYLTREE_FMT = "/users/ldog/alouis/ws2/GENOMICUS_SVN/data{0}/PhylTree.Ensembl.{0}.conf"
NEW_DUP_SUFFIX = re.compile(r'\.[A-Za-z`]+$')
ANCGENE_START = 'ENSGT'
ANCGENE2SP_PATTERN = r'([A-Z][A-Za-z_.-]+)(%s.*)$'
ANCGENE2SP = re.compile(ANCGENE2SP_PATTERN % ANCGENE_START)


#SPLIT_SPECIES_GENE = re.compile()

def print_if_verbose(*args, **kwargs):
    print(*args, **kwargs)


#def split_species_gene(nodename):
#    """Split a node name into its two parts (taxon + ancestral gene name)."""
#    try:
#        idx = nodename.index('ENS')
#    except ValueError:
#        try: # any subtring below here doesn't happen anymore with the last version.
#            idx = nodename.index('FBgn') # Drosophila
#        except ValueError:
#            try:
#                idx = nodename.index('WBGene') # Caenorhabditis
#            except ValueError:
#                try:
#                    idx = nodename.index('Y')
#                except ValueError:
#                    try:
#                        idx = nodename.index('Q0')
#                    except ValueError as err:
#                        err.args += ("ERROR: Invalid nodename %r" % nodename,)
#                        raise
#    return nodename[:idx].replace('.', ' '), nodename[idx:]
ucsc_conv_filename = '~/ws2/UCSC_genome_releases_full.tsv'
try:
    UCSC_CONVERSION = load_conversion(ucsc_conv_filename)
except FileNotFoundError:
    print("WARNING: conversion file not found: %r" % ucsc_conv_filename,
            file=sys.stderr)
    UCSC_CONVERSION = {}


def ultimate_seq2sp(seqname, ensembl_version=ENSEMBL_VERSION):
    """From a sequence name, find the corresponding species.
    Recognizes Ensembl gene IDs, Ensembl protein IDs, and also UCSC assembly
    names such as 'loxAfr3'"""
    try:
        sp = convert_gene2species(seqname, ensembl_version)
    except RuntimeError:
        try:
            sp = convert_prot2species(seqname, ensembl_version)
        except KeyError:
            assembly = re.match('[A-Za-z0-9]+', seqname).group()
            sp = UCSC_CONVERSION[assembly]
    return sp


def split_species_gene(nodename, ancgene2sp):
    match = ancgene2sp.match(nodename)
    try:
        taxon, genename = match.groups()
        taxon = taxon.replace('.', ' ')
    except AttributeError:
        taxon = None
        genename = nodename
    return taxon, genename


def name_missing_spe(parent_sp, ancestor, genename, parent_genename,
                     diclinks, ages=None): #, event='spe'):
    """given two taxa and the child gene name, return each missing speciation
    node inbetween."""
    try:
        ancestor_lineage = diclinks[parent_sp][ancestor]
    except KeyError as err:
        err.args = (err.args[0] +
                    "\nchild node : %s (%r)" % (genename, ancestor) +
                    "\nparent node: %s (%r)" % (parent_genename, parent_sp),)
        raise err
    # Links from diclinks must be corrected: unnecessary nodes removed, i.e
    # nodes with a single child and age 0. (introducing potential errors)
    # BUT: do not remove species (age = 0) !
    if ages:
        age_parent, age_child = ages[ancestor_lineage[0]], ages[ancestor_lineage[-1]]
        links = [link for link in ancestor_lineage[1:-1]
                     if age_parent >= ages[link] >= age_child and ages[link] > 0]
        ancestor_lineage = [ancestor_lineage[0]] + links + [ancestor_lineage[-1]]
    else:
        links = ancestor_lineage[1:-1]

    # If doesn't match species tree
    # same ancestor as child is possible for duplication node
    # So check for length 1 or 2
    new_node_names = []
    if links:
        # Alright to use parent_genename if not a duplication.
        # TODO: do not make the new name here, but in the `insert_...` function
        new_node_names = [link.replace(' ', '.') + parent_genename for link in links]

    ### TODO: refactor the part above to compute and check 'links' only once, 
    ###       during the following loop, and to check for *each* link age 
    ###       relatively to the previous one.

    if ages is not None:
        new_branch_dists = []
        total_len = ages[ancestor_lineage[0]] - ages[ancestor_lineage[-1]]
        for link_parent, link in zip(ancestor_lineage[:-1],
                                     ancestor_lineage[ 1:]):
            new_dist = float(ages[link_parent] - ages[link])
            if new_dist < 0:
                print("INVALID AGE:\n child node : %s (%r)\n" % (genename, ancestor),
                      "parent node: %s (%r)\n" % (parent_genename, parent_sp),
                      "total_len: %s - %s = %s\n" % (age_parent, age_child, total_len),
                      "new link: %s - %s\n" % (link_parent, link),
                      "new link ages: %s - %s" % (ages[link_parent], ages[link]),
                      file=sys.stderr)
                raise AssertionError('Tree %r: Invalid Age: %s (%s)' % \
                            (parent_genename, link_parent, ages[link_parent]))
            try:
                new_branch_dists.append(new_dist / total_len)
            except ZeroDivisionError as err:
                err.args += (parent_sp, ancestor, "please check ages.")
                print("WARNING:", err, file=sys.stderr)
                new_branch_dists.append(0)
    else:
        total_len = len(ancestor_lineage) - 1
        new_branch_dists = [1./total_len] * total_len

    return new_node_names, new_branch_dists, links
    # TODO: compute intermediate distances proportionally to the age of each
    # intermediate taxon. Will however be incorrect if one node is a duplication.


def insert_nodes(new_node_names, parent, child, new_taxa, new_dist_ratios=None):
    # conserve original distance
    if not new_dist_ratios:
        # split equally
        n_new_branches = len(new_node_names) + 1
        new_branch_dists = [float(child.dist) / n_new_branches] * n_new_branches
    else:
        new_branch_dists = [float(child.dist) * r for r in new_dist_ratios]

    new_node = parent
    new_nodes = []
    print_if_verbose("Inserted nodes\n (   %s\n  -> %s):" % (new_node.name, child.name))
    for new_name, new_taxon, new_dist in zip(new_node_names, new_taxa, new_branch_dists[:-1]):
        new_node = new_node.add_child(name=new_name, dist=new_dist)
        new_node.add_features(reinserted=True, S=new_taxon) #D="N"
        new_nodes.append(new_node)
        print_if_verbose(" -", new_name)
    #print_if_verbose()
    # Move the child on top of the created intermediate links
    new_node = new_node.add_child(child=child.detach(), dist=new_branch_dists[-1])
    return new_nodes


### TODO: write `insert_missing_dup`
def insert_missing_spe(parent_sp, ancestor, genename, parent_genename,
                       parent_node, child, diclinks, ages=None, event='spe'):
    """Insert missing speciation nodes between parent_node and child.
    
    If event='spe': parent_node is a speciation, 
    else if event='dup', it is a duplication and an extra speciation node
    must be appended."""
    new_node_names, new_branch_dists, new_taxa = name_missing_spe(parent_sp,
                                                        ancestor,
                                                        genename,
                                                        parent_genename,
                                                        diclinks,
                                                        ages)
    if new_node_names:
        print_if_verbose("Nodes to insert: %s" % new_node_names)
        return insert_nodes(new_node_names, parent_node, child, new_taxa, new_branch_dists)
        #return True
    # return false if no node was added
    #return False
    return []


def add_species_nodes_back(tree, diclinks, ages=None):
    """WRONG. DO NOT USE. Use `insert_species_nodes_back`
    Add missing species ancestors in gene tree"""
    # TODO: conserve branch length
    # Iterate from leaves to root
    for node in tree.iter_descendants("postorder"):
        if node.is_leaf():
            ancestor = convert_gene2species(node.name)
            genename = node.name
        else:
            ancestor, genename = split_species_gene(node.name, ancgene2sp)
        
        parent_node = node.up
        parent_ancestor, parent_genename = split_species_gene(parent_node.name, ancgene2sp)

        insert_missing_spe(parent_ancestor, ancestor, genename, parent_genename,
                           parent_node, node, diclinks, ages)


def suffixes_ok(parent, child, event):
    """parent: genename
        child: genename
        event: 'dup' or 'spe'"""
    #TODO: what to do with leaves (modern gene names)
    if event == 'spe':
        return parent == child
    elif event =='dup':
        new_suffix = child[len(parent):]
        return child.startswith(parent) and NEW_DUP_SUFFIX.match(new_suffix)
    else:
        raise RuntimeError("Invalid argument 'event' (must be 'dup' or 'spe')")


def suffix_count(parent, child):
    """count how many duplication suffixes were added between parent and child
    gene names. (suffixes like '.a', '.a.a', '.a.b'...)"""
    if not child.startswith(parent):
        raise RuntimeError("parent %r and child %r are not in the same lineage")

    difference = child[len(parent):]
    count = 0
    while NEW_DUP_SUFFIX.search(difference):
        difference = NEW_DUP_SUFFIX.sub('', difference)
        count += 1
    return count


def suffix_list(parent, child):
    """list duplication suffixes that were added between parent and child
    gene names. (suffixes like '.a', '.b', '.`b', '.ag'...)"""
    if not child.startswith(parent):
        #raise RuntimeError("parent %r and child %r are not in the same lineage")
        return None

    difference = child[len(parent):]
    suffixes = []
    match = NEW_DUP_SUFFIX.search(difference)
    while match:
        suffixes.append(match.group())
        difference = NEW_DUP_SUFFIX.sub('', difference)
        match = NEW_DUP_SUFFIX.search(difference)
    return suffixes


def get_mrca(parent_sp, children_sp, diclinks):
    """Get most recent common ancestor of all children species, given a root
    'parent_sp'."""
    children_anc = [diclinks[parent_sp][ch_sp] for ch_sp in children_sp]
    for next_parents in zip(*children_anc):
        #print(parent_sp, next_parents)
        if len(set(next_parents)) > 1:
            break
        else:
            mrca = next_parents[0]
    try:
        return mrca
    except UnboundLocalError as err:
        err.args = (err.args[0] + " (%s, %s)" %(parent_sp, children_sp), )
        raise


def insert_species_nodes_back(tree, ancgene2sp, diclinks, ages=None,
                              fix_suffix=True, force_mrca=False,
                              ensembl_version=ENSEMBL_VERSION, treebest=False):
    if treebest:
        print_if_verbose("  Reading from TreeBest format")
        get_species    = lambda node: (node.S.replace('.', ' '), node.name.split('_')[0])
        split_ancestor = lambda node: (node.S.replace('.', ' '), node.name) 
    else:
        get_species = lambda node: (ultimate_seq2sp(node.name, ensembl_version),
                                    node.name)
        split_ancestor = lambda node: split_species_gene(node.name, ancgene2sp)

    print_if_verbose("* Insert missing nodes:")
    ### Insert childs *while* iterating.
    ### Beware not to iterate on these new nodes:
    ### iterate from leaves to root ('postorder') and insert between node and child
    for node in tree.traverse('postorder'):
        print_if_verbose(" * %r" % node.name)
        # Copying this list is CRUCIAL: children get inserted in node.children,
        # so you will iterate over your inserted children otherwise.
        node_children = copy(node.children)
        if node_children:
            parent_sp, parent_gn = split_ancestor(node)
            ### Delete this node if species is not recognized
            if parent_sp is None:
                print("WARNING: taxon not found in (%r, %r). Deleting node." %
                        (parent_sp, parent_gn), file=sys.stderr)
                node.delete(prevent_nondicotomic=False,
                            preserve_branch_length=True)
                continue

            child_sp = []
            child_gn = []
            nodes_without_taxon = []
            for child in node_children:
                print_if_verbose("  - child %r" % child.name)
                if child.is_leaf():
                    try:
                        ancestor, genename = get_species(child)
                    except KeyError as err:
                        print("WARNING: Leaf %r not in an extant species" % \
                              child.name,
                              file=sys.stderr)
                        ancestor, genename = split_ancestor(child)
                else:
                    ancestor, genename = split_ancestor(child)
                if ancestor not in diclinks:
                    print("WARNING: taxon %r absent from phylogeny. "
                          "Deleting node %r." % (ancestor, genename),
                          file=sys.stderr)
                    nodes_without_taxon.append(child)
                    child.delete(prevent_nondicotomic=False,
                                 preserve_branch_length=True)
                else:#    continue
                    child_sp.append(ancestor)
                    child_gn.append(genename)

            
            # MUST UPDATE because the children list just changed.
            for n_without_taxon in nodes_without_taxon:
                node_children.remove(n_without_taxon)

            # Tricky part: possible situations:
            # tree : parent X -> children Y1, Y2 (Y3...)
            # 1. All equal     -> duplication
            # 2. All different -> speciation (but must check that X is the MRCA)
            # 
            # 3. At least one child = parent, but not all
            #   -> a duplication (+ missing speciation)
            # 
            # 4. At least one pair of child is equal, but they are all
            #    different from parent.
            #   -> speciations after X are missing (should not happen)
            #   + a duplication is missing.
            #
            # It seems like cases 2 and 4 can be grouped together (spe) as 
            # "at least one pair of children different" -> check MRCA fits
            #
            # 1 and 3 can be grouped as duplication: ""

            # Compare the parent with its descendants
            compare_vertical = [parent_sp == sp for sp in child_sp]
            # compare each pair of children (True if at least one pair)
            compare_horizontal = len(child_sp) > len(set(child_sp))

            # 1. all(compare_vertical) and all(compare_horiz)
            # 2. not any(compare_vertical) and not any(compare_horiz)
            # 3. any(compare_vertical) and not all(compare_vertical)
            # 4. any(compare_horizontal) and not any(compare_vertical)

            #if all(test_same_species):
            #    # duplication.
            #    # check suffixes
            #    event = 'dup'
            #elif not any(test_same_species):
            #    # speciation
            #    # check suffixes
            #    event = 'spe'
            if any(compare_vertical):
                event = 'dup'
            else:
                event = 'spe'
                # Check that parent is the MRCA
                mrca = get_mrca(parent_sp, child_sp, diclinks)
                if not parent_sp == mrca:
                    # I didn't expect this to happen. Let's fix it (not raise)
                    print("WARNING: Unexpected case: parent %r is not the "
                          "MRCA of %s." % (parent_sp, child_sp),
                          end=' ', file=sys.stderr)
                    if not force_mrca:
                        print("Ignoring.", file=sys.stderr)
                        # Must change event to dup.
                        event = 'dup'
                    else:
                        print("Setting common parent to MRCA %r." % mrca,
                              file=sys.stderr)
                        # Choose the length of the common branch to add according
                        # to the closest child (shortest branch)
                        closest_child = min(node_children, key=lambda node: node.dist)
                        closest_i = node_children.index(closest_child)
                        print_if_verbose('  --- closest child %r:' % closest_child.name)
                        new_nodes = insert_missing_spe(parent_sp,
                                            child_sp[closest_i],
                                            child_gn[closest_i],
                                            parent_gn,
                                            node, closest_child,
                                            diclinks, ages)
                        mrca_i, mrca_node = [(i, n) for i,n in enumerate(new_nodes)
                                                if n.name.startswith(
                                                        mrca.replace(' ', '.'))][0]
                        print_if_verbose('MRCA node: %r' % mrca_node)
                        mrca_dist = sum(n.dist for n in new_nodes[1:(mrca_i+1)])
                        for i, other_child in enumerate(node_children):
                            print_if_verbose('  --- other child %r:' % other_child.name,
                                             end=' ')
                            if other_child != closest_child:
                                print_if_verbose('re-branch to the MRCA: %s ->'
                                                    % other_child.up.name,
                                                 end=' ')
                                mrca_node.add_child(child=other_child.detach(),
                                                    dist=other_child.dist - mrca_dist)
                                print_if_verbose(other_child.up.name)
                                #insert_missing_spe(mrca, child_sp[i], child_gn[i],
                                #                   parent_gn,
                                #                   mrca_node, other_child,
                                #                   diclinks, ages)
                            else:
                                print_if_verbose()
                        print_if_verbose(node.get_ascii())

                        node_children.remove(closest_child)
                        child_sp = child_sp[:closest_i] + child_sp[(closest_i+1):]
                        child_gn = child_gn[:closest_i] + child_gn[(closest_i+1):]
                        node = mrca_node
                        parent_sp = mrca

                if compare_horizontal:
                    #if all(compare_horizontal):
                    #    raise AssertionError("Descendants are duplicates, but "
                    #                         "the parent node species differs.")
                    #    # Same as the previous assertion (i.e. is not MRCA)
                    #else:
                    raise AssertionError("Unexpected case: some but not "
                                         "all descendants are duplicates.")

            if event == 'dup':
                if not all(compare_vertical):
                    print_if_verbose("implicit dup+spe")
                    # there is a dup + a spe. Need to add nodes
                    ### TODO: check this missing node **before** this node:
                    ###       if node_taxon != parent_taxon and \
                    ###           sisternode_taxon == parent_taxon
                    ###       HOWEVER this will potentially mess the new_dists
                    added_suffixes = [suffix_list(parent_gn, gn) for gn in child_gn]
                    added_suffixes = [suf[0].lstrip('.') if suf else None for suf in added_suffixes]
                    
                    ### TODO: put this section into a `fix_suffix` function.
                    suff_count = 1 #len(set((suf for suf in added_suffixes if suf)))
                    for i, (child, ancestor, genename, suf) in \
                            enumerate(zip(node_children, child_sp, child_gn,
                                            added_suffixes)):
                        if ancestor != parent_sp:
                            # append proper suffix if there isn't one.
                            if not suf and fix_suffix:
                                suf = chr(96+suff_count)
                                while suf in added_suffixes:
                                    suff_count += 1
                                    suf = chr(96+suff_count)
                                added_suffixes.append(suf)
                                genename = parent_gn + '.' + suf
                                child_gn[i] = genename
                            # Insert back the needed speciation event after.
                            # (make this node a duplication again)
                            spe_node_name = parent_sp.replace(' ', '.') + genename
                            print_if_verbose("Inserted node (spe after dup):",
                                             spe_node_name)
                            # Then check the succession of ancestors
                            new_node_names, new_branch_dists, new_taxa = \
                                    name_missing_spe(parent_sp,
                                                     ancestor,
                                                     genename,
                                                     genename, #because speciation
                                                     diclinks,
                                                     ages)
                            new_node_names.insert(0, spe_node_name)
                            new_taxa.insert(0, parent_sp)
                            # adjust the length of this dup-spe branch
                            n_new_nodes = len(new_node_names)
                            new_branch_dists = [1./(n_new_nodes+1)] + \
                                               [r * n_new_nodes/(n_new_nodes+1)
                                                   for r in new_branch_dists]
                            new_nodes = insert_nodes(new_node_names, node,
                                                     child, new_taxa, new_branch_dists)
                            # MUST update node_children after the previous edit
                            node_children[i] = new_nodes[0]
                            child_sp[i] = parent_sp
            
            #node_children = copy(node.children)
            for child, ancestor, genename in zip(node_children, child_sp, child_gn):
                if not child.is_leaf():
                    try:
                        assert suffixes_ok(parent_gn, genename, event) or not fix_suffix
                    except AssertionError:
                        print("WARNING: inconsistent suffixes:\n",
                              "  - node            %r\n" % node.name,
                              "  - children        %s\n" % (node_children,),
                              "  - parent sp       %r\n" % parent_sp,
                              "  - species         %r\n" % ancestor,
                              "  - parent genename %r\n" % parent_gn,
                              "  - genename        %r\n" % genename,
                              "  - event           %r" % event, file=sys.stderr)
                        #raise

                # Nothing gets inserted if it is a duplication.
                insert_missing_spe(parent_sp, ancestor, genename, parent_gn,
                                   node, child, diclinks, ages, event=event)


def search_by_ancestorlist(tree, ancestorlist, latest_ancestor=False):
    """Return an iterator over all basal nodes belonging to an ancestor in 
    ancestorlist. If `latest_ancestor` is True, return the most recent nodes
    belonging to one of these ancestors."""
    if not latest_ancestor:
        #print("")
        def stop_at_any_ancestor(node):
            return any(node.name.startswith(anc) for anc in ancestorlist)
    else:
        def stop_at_any_ancestor(node):
            match_anc = [anc for anc in ancestorlist if node.name.startswith(anc)]
            if match_anc:
                anc = match_anc[0]
                # test for any
                return node.is_leaf() or \
                       not all(ch.name.startswith(anc) for ch in node.children)
            else:
                return False
    
    return tree.iter_leaves(is_leaf_fn=stop_at_any_ancestor)


def with_dup(leafnames):
    #leafspecies = [convert_gene2species(leaf) for leaf in leafnames]
    leafspecies = [ultimate_seq2sp(leaf) for leaf in leafnames]
    return (len(leafspecies) > len(set(leafspecies)))


def get_basal(nodes, maxsize):
    """Identify `maxsize` most basal nodes from a list of sister nodes.

    Return 2 lists:
        - selected basal nodes,
        - excluded nodes (to detach)
    
    Basal means: closer to the root (topology only, i.e in number of nodes).
    
    - `nodes`: a *list* of TreeNode instances;
    - `maxsize`: integer.

    Return ([],[]) if maxsize >= number of leaves.
    """
    
    get_dist = lambda n: n.dist
    nodes.sort(key=get_dist, reverse=True)

    # Not enough leaves
    if sum(len(n) for n in nodes) <= maxsize:
        return [], []
        #return list(chain(n.get_leaves() for n in nodes))

    # Traverse in a `levelorder` strategy (breadth-first-search)
    # until the limit size is reached.
    
    #kept = []
    #while len(kept) + len(nodes) < maxsize:
    #    # Remove closest node and add its children. Keep if no children.
    #    node = nodes.pop(0)
    #    if node.children:
    #        nodes.extend(node.children)
    #    else:
    #        kept.append(node)

    #base = []
    #while True:
    #    for node in nodes:
    #        base.extend(node.children or [node])
    #        if len(base) >= maxsize:
    #            return base[:maxsize]
    #    # start at next level
    #    nodes = base
    #    base = []

    kept = 0
    while len(nodes) < maxsize:
        nextnodes = nodes[kept].children
        if nextnodes:
            nodes.pop(kept)
            nodes.extend(sorted(nextnodes, key=get_dist, reverse=True))
        else:
            kept += 1

    #return (to keep, to detach)
    return nodes[:maxsize], nodes[maxsize:]


def reroot_with_outgroup(node, maxsize=0):
    """Goes up the tree (towards the root) until it finds outgroup taxa.
    
    - Only keep at most `maxsize` leaves in the outgroup.
    - Keep all leaves if maxsize < 0.
    - Return None if there is no outgroup.

    If the common ancestor with the outgroup is a multifurcation,
    process all of the sister taxa as a single outgroup.

    Select outgroup leaves to keep that are the most distant to each other
    (having the most basal latest common ancestor).
    """
    if maxsize == 0:
        return node
    
    root = node
    while len(root) == len(node):
        if root.is_root():
            print("WARNING: no outgroup available for node %r" % node.name,
                  file=sys.stderr)
            
            return
        # Else go to parent node
        node = root
        root = node.up

    # If needed, reduce the size of the outgroup
    # clade to the specified number.
    outgroups = node.get_sisters()

    if maxsize > 0 and \
            sum(len(outgroup) for outgroup in outgroups) > maxsize:

        for outgroup in outgroups:
            #outgroup.ladderize()
            # Now the first listed leaves are the closest.

            # Deepcopy slows it down too much. Alternate solution:
            # Add a uniq mark, because I'm not sure that searching on name
            # would be unambiguous.
            outgroup.add_feature('tmpmark', True)

        # MAKE A COPY (because this part could be reused later)
        # Need to use `deepcopy` to preserve the ancestor information.
        #outgroup = deepcopy(outgroup)
        #root = outgroup.up
        root = root.copy()

        base, todetach = get_basal(root.search_nodes(tmpmark=True), maxsize)
        
        # Now that we know the basal tree:
        # keep only one leaf per identified basal node.
        for basal in base:
            basal.prune([basal.get_closest_leaf()[0]], preserve_branch_length=True)
            # WARN: 'prune' deletes all intermediate nodes leading to the kept leaves
        
        # And detach the extraneous basal nodes
        for nodetodetach in todetach:
            nodetodetach.detach()

        for outgroup in outgroups:
            outgroup.del_feature('tmpmark')

    return root


def save_subtrees(treefile, ancestorlists, ancestor_regexes, ancgene2sp,
        diclinks, treebest=False, ages=None, fix_suffix=True, force_mrca=False,
        latest_ancestor=False, ensembl_version=ENSEMBL_VERSION, outdir='.',
        only_dup=False, one_leaf=False, outgroups=0, dry_run=False):
    #print_if_verbose("* treefile: " + treefile)
    #print("treebest = %s" % treebest, file=sys.stderr)
    outfiles_set = set() # check whether I write twice to the same outfile
    if treefile == '-': treefile = '/dev/stdin'
    try:
        tree = ete3.Tree(treefile, format=1)
    except ete3.parser.newick.NewickError as err:
        err.args = (err.args[0] + 'ERROR with treefile %r' % treefile,)
        raise
    output_features = tree.features.union(('reinserted',)) - set(('name', 'support', 'dist'))
    insert_species_nodes_back(tree, ancgene2sp, diclinks, ages, fix_suffix,
                              force_mrca, ensembl_version, treebest)
    print_if_verbose("* Searching for ancestors:")
    for ancestor, ancestorlist in ancestorlists.items():
        print_if_verbose(ancestor)
        ancestor_regex = ancestor_regexes[ancestor]
        for ancestornodeid, node in enumerate(search_by_ancestorlist(tree, ancestorlist, 
                                                        latest_ancestor)):
        #for ancestornodeid, node, root in enumerate(search_by_ancestorlist(tree, ancestorlist, 
        #                                                    latest_ancestor)):
            leafnames = node.get_leaf_names()
            #print(node.name)
            #print(node.get_ascii())
            if len(leafnames) > 1 or one_leaf:
                # check that there is at least one duplication
                if not only_dup or with_dup(leafnames):
                    ### TODO: when you *know* there can be duplicated node
                    ###       names, use a custom unique id for each node.
                    outname = ancestor_regex.sub(ancestor, node.name)
                    if outname == ancestor:
                        outname += "%02d" % ancestornodeid
                    #outname = re.sub('^[A-Za-z_.-]+(?=ENSGT)', node.name, ancestor)
                    #print_if_verbose(('generate outname:\n'
                    #                  'regex: %r\n'
                    #                  'node:  %r\n'
                    #                  'anc:   %r\n'
                    #                  'outname: %r') % (ancestor_regex.pattern,
                    #                                  node.name,
                    #                                  ancestor,
                    #                                  outname))

                    outfile = os.path.join(outdir, outname + '.nwk')
                    if outdir == '-':
                        outfile = None
                    elif outfile in outfiles_set:
                        # Not sure this case can happen, but better prevent it.
                        # Months later: this case happened.
                        print("ERROR: Cannot output twice to %r" % outfile,
                              file=sys.stderr)
                        sys.exit(1)
                    
                    # Add requested outgroups:
                    root = reroot_with_outgroup(node, outgroups)
                    if not root:
                        continue

                    if dry_run:
                        print("node %r (%s)\n\\-> %s" % (node.name, ancestor, outfile))
                    else:
                        print_if_verbose("Writing to %r." % outfile)
                        outtree = root.write(format=1,
                                             format_root_node=True,
                                             outfile=outfile,
                                             features=output_features)
                        if outtree is not None:
                            print(outtree)
                    outfiles_set.add(outfile)
    #print_if_verbose("Output %d trees." % len(outfiles_set))
    return list(outfiles_set)


def save_subtrees_process(params):
    print("* Input tree: %r" % params[0], file=sys.stderr)
    ignore_errors = params.pop()
    try:
        outfiles = save_subtrees(*params)
    except BaseException as err:
        if ignore_errors:
            print("Ignore %r: %r" % (params[0], err), file=sys.stderr)
            return 0, []
        else:
            raise
    return 1, outfiles # return a value to let Pool.map count the number of results.


def parallel_save_subtrees(treefiles, ancestors, ncores=1, outdir='.',
                           outsub=None, treebest=False, only_dup=False,
                           one_leaf=False, fix_suffix=True, force_mrca=False,
                           latest_ancestor=False, outgroups=0, dry_run=False,
                           ignore_errors=False,
                           ensembl_version=ENSEMBL_VERSION):
    ### WARNING: uses global variables here, that are changed by command line
    phyltree = PhylTree.PhylogeneticTree(PHYLTREE_FMT.format(ensembl_version))
    # Crucial point: the pattern alternatives must be sorted by age, so that
    # you don't match an ancestor whose name is contained in its descendant name
    # (like theria is contained in eutheria)
    ancgene2sp = re.compile(r'('
                            + r'|'.join(list(phyltree.listSpecies) + 
                                        sorted(phyltree.listAncestr,
                                               key=lambda a:len(a),
                                               reverse=True)).replace(' ','\.')
                            + r')(.*)$')
    ancestorlists = {}
    ancestor_regexes = {}
    for anc_lowercase in ancestors:
        ancestor = anc_lowercase.capitalize()
        ancestorlist = sorted(phyltree.allDescendants[ancestor],
                              key=lambda anc: -phyltree.ages[anc])
        ancestorlists[anc_lowercase] = ancestorlist
        ancestor_regexes[anc_lowercase] = re.compile('^(%s)' % \
                                    '|'.join(ancestorlist).replace(' ', '.'))

    diclinks = phyltree.dicLinks.common_names_mapper_2_dict()
    ages = phyltree.ages.common_names_mapper_2_dict()

    if outsub:
        def format_outdir(treefile):
            return outdir.format(''.join(os.path.basename(treefile).split(outsub)[:-1]))

    else:
        def format_outdir(treefile):
            return outdir.format(os.path.splitext(os.path.basename(treefile))[0])

    generate_args = [[treefile,
                      ancestorlists,
                      ancestor_regexes,
                      ancgene2sp,
                      diclinks,
                      treebest,
                      ages,
                      fix_suffix,
                      force_mrca,
                      latest_ancestor,
                      ensembl_version,
                      format_outdir(treefile),
                      only_dup,
                      one_leaf,
                      outgroups,
                      dry_run,
                      ignore_errors] for treefile in treefiles]

    n_input = len(treefiles)
    print("To process: %d input trees" % n_input, file=sys.stderr)
    if ncores > 1:
        pool = mp.Pool(ncores)
        outputs = pool.map(save_subtrees_process, generate_args)
        return_values, outfiles = zip(*outputs)
        progress = sum(return_values)
        outfiles = [outf for outf_group in outfiles for outf in outf_group]
    else:
        progress = 0
        outfiles = []
        for args in generate_args:
            return_value, some_outfiles = save_subtrees_process(args)
            progress += return_value
            outfiles.extend(some_outfiles)
    print("Finished processing %d/%d input trees. Output %d trees." % \
            (progress, n_input, len(outfiles)), file=sys.stderr)
    return outfiles


def parse_treefiles(treefiles_file):
    with open(treefiles_file) as stream:
        return [line.rstrip() for line in stream]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("treefile")
    parser.add_argument("ancestors", nargs='+')
    parser.add_argument("--fromfile", action="store_true", help="read treefile"\
                        " names from the first argument")
    parser.add_argument("-o", "--outdir", default='./{0}', help="[%(default)s]")
    parser.add_argument("-s", "--outsub", help="alternative splitting " \
                        "character to remove the extension from the basename "\
                        "of the treefile (used by '{0}' in --outdir).")
    parser.add_argument('-t', '--treebest', action='store_true',
                        help='input trees are in TreeBest output format ' \
                             '[Experimental!]')
    parser.add_argument("--only-dup", action="store_true",
                        help="do not extract trees that don't have at least "\
                             "one duplication")
    parser.add_argument("--one-leaf", action="store_true",
                        help="also output trees that have only one leaf.")
    parser.add_argument("-e", "--ensembl-version", type=int, default=85,
                        help="[%(default)s]")
    parser.add_argument("-p", "--phyltree-fmt", default=PHYLTREE_FMT,
                        help="Phylogenetic species tree "\
                        "in LibsDyogen PhylTree format. Can contain the string"\
                        " '{0}' which will be replaced by the Ensembl version"\
                        " [%(default)s]")
    #parser.add_argument("-a", "--ancgene-start", default="ENSGT", help="start"\
    #                    "of ancgenes (regex) [%(default)s]")
    parser.add_argument("--nofix-suffix", action="store_false",
                        dest="fix_suffix", help="disable attempt to fix gene "\
                        "suffixes (when adding back duplicated nodes)")
    parser.add_argument("--force-mrca", action="store_true",
                        help="If parent is not the species MRCA of its " \
                             "children, rebranch children on the MRCA.")
    parser.add_argument("-l", "--latest-ancestor", "--latest",
                        action="store_true",
                        help="When an ancestor has multiple paralogs of the " \
                        "gene, output one tree per paralog and root each tree"\
                        " at the latest gene, (instead of rooting at the " \
                        "ancestral gene)")
    parser.add_argument("--outgroups", metavar='S', type=int, default=0, 
                        help="Save the subtree including the outgroup sister "\
                             "clade of maximum size %(metavar)s. Set to '-1' "\
                             " to keep the full clade [Default S=%(default)s]."\
                             "NOTE: with `-l`, the outgroup can be a paralog "\
                             "tree, otherwise it's necessarily in sister species.")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="only print out the output files it would produce")
    parser.add_argument("--ncores", type=int, default=1, help="Number of cores")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-i", "--ignore-errors", action="store_true", 
                        help="On error, print the error and continue the loop.")
    #parser.add_argument("-O", "--to-stdout", action="store_true")
    
    args = parser.parse_args()
    dargs = vars(args)

    PHYLTREE_FMT = dargs.pop("phyltree_fmt")
    #ENSEMBL_VERSION = dargs.pop("ensembl_version")
    #ANCGENE_START = dargs.pop("ancgene_start")
    #ANCGENE2SP = re.compile(ANCGENE2SP_PATTERN % ANCGENE_START)

    if not dargs.pop("verbose"):
        # mute the custom print function
        def print_if_verbose(*args, **kwargs):
            pass

    if dargs.pop("fromfile"):
        treefiles = parse_treefiles(dargs.pop("treefile"))
    else:
        treefiles = [dargs.pop("treefile")]
    parallel_save_subtrees(treefiles, **dargs)

