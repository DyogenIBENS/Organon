#!/usr/bin/env python3

from __future__ import print_function

"""Extract all gene subtrees rooted at a given ancestral taxon.
Also add missing speciation nodes."""


import sys  # Changing stdout for all the interpreter only works when using sys.stdout
from sys import stdin, stderr, exit, setrecursionlimit
from io import StringIO
import os.path as op
import re
from functools import partial

try:
    import argparse_custom as argparse
    #print('ARGPARSE CUSTOM')
except ImportError:
    import argparse

import logging
import multiprocessing as mp
try:
    from multiprocessing_logging import install_mp_handler
except ImportError:
    logger.warning('multiprocessing_logging module not found, your stderr logs will likely be messed up.')
    def install_mp_handler(logger):
        pass

from copy import copy

import ete3
import LibsDyogen.myPhylTree as PhylTree

from dendro.parsers import read_multinewick, iter_from_ete3
from genomicustools.identify import ultimate_seq2sp
from dendro.bates import iter_distleaves
from dendro.trimmer import thin_ete3 as thin

logger = logging.getLogger(__name__)

#stdoutlog = logging.getLogger(__name__ + '.stdout')
#stdouth = logging.StreamHandler(stdout)
#stdouth.setFormatter(logging.Formatter("%(message)s"))
#stdoutlog.addHandler(stdouth)
#stdoutlog.setLevel(logging.INFO)


ENSEMBL_VERSION = 85
PHYLTREE_FMT = "/users/ldog/alouis/ws2/GENOMICUS_SVN/data{0}/PhylTree.Ensembl.{0}.conf"
NEW_DUP_SUFFIX = re.compile(r'\.([A-Z@]+|[a-z`@]+)$')  # The backtick is allowed because of a bug in older version of LibsDyogen.myProteinTree.getDupSuffix
NEW_ROOTDUP_SUFFIX = re.compile(r'\.[A-Z@]+$')
ANCGENE_START = 'ENSGT'
ANCGENE2SP_PATTERN = r'([A-Z][A-Za-z_.-]+)(%s.*)$'
ANCGENE2SP = re.compile(ANCGENE2SP_PATTERN % ANCGENE_START)
RENAME_ROOT_SUBST = r's/^({descendants})(?=ENSGT|\b)/{ancestor}/'  # used to replace '{descendants}' by `ancestor` in ancestor_regexes. Customise in order to avoid duplicated outnames.

#SPLIT_SPECIES_GENE = re.compile()

def parse_substitution_expr(subst):
    """Parse sed-like substitution expression: s/old/new/.
    
    You can also use another delimiter, for example '%': s%old%new%
    """
    assert subst.startswith('s')
    delim = subst[1]
    _, pattern, repl, _ = subst.split(delim) # fails (good thing) if wrong number of delimiters
    return pattern, repl

def print_if_verbose(*args, **kwargs):
    print(*args, **kwargs)
    #stdoutlog.info(*args, **kwargs)


def split_species_gene(nodename, ancgene2sp):
    match = ancgene2sp.match(nodename)
    try:
        taxon, genename = match.groups()
        taxon = taxon.replace('.', ' ')
    except AttributeError:
        taxon = None
        genename = nodename
    return taxon, genename


def parse_species_genename(child, get_species, split_ancestor):
    if child.is_leaf():
        try:
            ancestor, genename = get_species(child)
        except KeyError as err:
            logger.warning("Leaf %r not in an *extant* species",
                           child.name)
            ancestor, genename = split_ancestor(child)
    else:
        ancestor, genename = split_ancestor(child)
    return ancestor, genename


def name_missing_spe(parent_sp, ancestor, genename, parent_genename,
                     diclinks, ages=None): #, event='spe'):
    # ~~> dendro.reconciled
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
                logger.error(
                     "\nINVALID AGE:\n child node : %s (%r)\n" % (genename, ancestor)
                     + "parent node: %s (%r)\n" % (parent_genename, parent_sp)
                     + "total_len: %s - %s = %s\n" % (age_parent, age_child, total_len)
                     + "new link: %s - %s\n" % (link_parent, link)
                     + "new link ages: %s - %s" % (ages[link_parent], ages[link]))
                raise AssertionError('Tree %r: Invalid Age: %s (%s)' % \
                            (parent_genename, link_parent, ages[link_parent]))
            try:
                new_branch_dists.append(new_dist / total_len)
            except ZeroDivisionError as err:
                err.args += (parent_sp, ancestor, "please check ages.")
                if link_parent != link:
                    logger.warning(err)
                else:
                    #assert
                    pass
                new_branch_dists.append(0)
    else:
        total_len = len(ancestor_lineage) - 1
        new_branch_dists = [1./total_len] * total_len

    return new_node_names, new_branch_dists, links
    # TODO: compute intermediate distances proportionally to the age of each
    # intermediate taxon. Will however be incorrect if one node is a duplication.


def insert_nodes(new_node_names, parent, child, new_taxa, new_dist_ratios=None):
    # ~~> dendro.reconciled
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
        if getattr(child, 'P', None) is not None:
            new_node.add_feature('P', child.P)
        new_nodes.append(new_node)
        print_if_verbose(" - %s" % new_name)
    #print_if_verbose()
    # Move the child on top of the created intermediate links
    new_node = new_node.add_child(child=child.detach(), dist=new_branch_dists[-1])
    # new_nodes.append(new_node)?
    return new_nodes


### TODO: write `insert_missing_dup`
def insert_missing_spe(parent_sp, ancestor, genename, parent_genename,
                       parent_node, child, diclinks, ages=None, event='spe'):
    # ~~> dendro.reconciled
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


def suffixes_ok(parent, child, event): # ~~> genomicustools/dendro?
    """parent: *genename* (Better if does not contain the taxon name);
        child: *genename*
        event: 'dup' or 'spe'
    """
    #TODO: what to do with leaves (modern gene names)
    if event == 'spe':
        return parent == child
    elif event =='dup':
        expected_suffix = NEW_DUP_SUFFIX
        if not child.startswith(parent):
            # Might be a duplication at the root of the tree, in which case
            # the suffix should be upper-case letters.
            assert parent.find(ANCGENE_START) == 0 and child.find(ANCGENE_START) == 0
            parent = parent.split('.')[0]  # Keep only the 'tree_name'.
            expected_suffix = NEW_ROOTDUP_SUFFIX
            # TOCHECK: that the taxon is the root taxon for this proteintree, 
            #          and that the child node is a speciation.

        if not child.startswith(parent):
            return False
            
        new_suffix = child[len(parent):]
        return expected_suffix.match(new_suffix)

    else:
        raise ValueError("Invalid argument 'event' (must be 'dup' or 'spe')")


def suffix_count(parent, child):
    """count how many duplication suffixes were added between parent and child
    gene names. (suffixes like '.a', '.a.a', '.a.b'...)"""
    if not child.startswith(parent):
        raise ValueError("parent %r and child %r are not in the same lineage")

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
        #raise ValueError("parent %r and child %r are not in the same lineage")
        return None

    difference = child[len(parent):]
    suffixes = []
    match = NEW_DUP_SUFFIX.search(difference)
    while match:
        suffixes.append(match.group())
        difference = NEW_DUP_SUFFIX.sub('', difference)
        match = NEW_DUP_SUFFIX.search(difference)
    return suffixes


def get_mrca(parent_sp, children_sp, diclinks): # ~~> a myPhylTree annex?.
    # Already in dicParents
    """Get most recent common ancestor of all children species, given a root
    'parent_sp'."""
    children_anc = [diclinks[parent_sp][ch_sp] for ch_sp in children_sp]
    for next_parents in zip(*children_anc):  #FIXME NOT OK
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
    #return phyltree.lastCommonAncestor(children_sp)


def insert_species_nodes_back(tree, parse_species_genename, diclinks, ages=None,
                              fix_suffix=True, force_mrca=False):

    print_if_verbose("* Insert missing nodes:")
    ### Insert childs *while* iterating.
    ### Beware not to iterate on these new nodes:
    ### iterate from leaves to root ('postorder') and insert between node and child
    for node in tree.traverse('postorder'):
        print_if_verbose(" * %r" % node.name)
        # Copying this list is CRUCIAL: children get inserted in node.children,
        # so you will iterate over your inserted children otherwise.
        node_children = node.get_children()
        if node_children:
            parent_sp, parent_gn = parse_species_genename(node)
            ### Delete this node if species is not recognized
            if parent_sp is None:
                logger.warning("Taxon not found in (%r, %r). Deleting node.",
                               parent_sp, parent_gn)
                for child in node.children:
                    child.dist += node.dist
                node.delete(prevent_nondicotomic=False,
                            preserve_branch_length=False)

                continue

            child_sp = []
            child_gn = []
            nodes_without_taxon = []
            for child in node_children:
                print_if_verbose("  - child %r" % child.name)
                ancestor, genename = parse_species_genename(child)
                if ancestor not in diclinks:
                    logger.warning("Taxon %r absent from phylogeny. "
                                   "Deleting node %r.", ancestor, genename)
                    nodes_without_taxon.append(child)
                    for nextchild in child.children:
                        nextchild.dist += child.dist
                    child.delete(prevent_nondicotomic=False,
                                 preserve_branch_length=False)
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
            same_taxa_vertically = [parent_sp == sp for sp in child_sp]
            # compare each pair of children (True if at least one pair)
            any_same_taxa_horizontally = len(child_sp) > len(set(child_sp))

            # 1. all(same_taxa_vertically) and all(compare_horiz)
            # 2. not any(same_taxa_vertically) and not any(compare_horiz)
            # 3. any(same_taxa_vertically) and not all(same_taxa_vertically)
            # 4. any(same_taxa_horizontally) and not any(same_taxa_vertically)

            #if all(test_same_species):
            #    # duplication.
            #    # check suffixes
            #    event = 'dup'
            #elif not any(test_same_species):
            #    # speciation
            #    # check suffixes
            #    event = 'spe'
            #event = infer_gene_event(node, parent_sp, child_sp)
            if any(same_taxa_vertically):
                event = 'dup'
            else:
                event = 'spe'
                node.add_feature('event', 'spe')
                # Check that parent is the MRCA
                if len(child_sp)>1:
                    mrca = get_mrca(parent_sp, child_sp, diclinks)
                else:
                    mrca = parent_sp
                if parent_sp != mrca:
                    # I didn't expect this to happen. Let's fix it (not raise)
                    msg = ("Unexpected case: parent %r is not the MRCA of %s. "
                            % (parent_sp, child_sp))
                    node.add_feature('event', 'spe not mrca')
                    if any_same_taxa_horizontally:
                        # Must change event to dup.
                        event = 'dup'
                        node.add_feature('D', 1)  # True or "dubious"?
                        node.add_feature('event', 'dup not mrca')

                    if not force_mrca:
                        logger.warning(msg + "Ignoring.")
                    else:
                        ### TODO: when several consecutive nodes have this
                        ###       problem, the fix must be applied to all
                        ###       (recursively?)
                        logger.warning(msg + "Setting common parent to MRCA %r.",
                                       mrca)
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

                # We detect a speciation. Double-check
                if any_same_taxa_horizontally:
                    msg = "Unexpected case: some but not all descendants are "\
                          "duplicates: %s: %s (%s)" % \
                                (node.name, node_children, event)
                    node.add_feature('event', 'spe+dup')
                    logger.error(msg)
                    #raise RuntimeError(msg)
                    #if all(same_taxa_horizontally):
                    #    raise AssertionError("Descendants are duplicates, but "
                    #                         "the parent node species differs.")
                    #    # Same as the previous assertion (i.e. is not MRCA)
                    #else:
                    #raise AssertionError()

            if event == 'dup':
                node.add_feature('event', "dup")
                if not all(same_taxa_vertically):
                    print_if_verbose("implicit dup+spe")
                    node.event = "implicit dup+spe"
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
                            ###TODO: do not raise the warning about clade ages

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
                        logger.warning("Inconsistent suffixes:\n"
                              + "  - node            %r\n" % node.name
                              + "  - children        %s\n" % (node_children,)
                              + "  - parent sp       %r\n" % parent_sp
                              + "  - species         %r\n" % ancestor
                              + "  - parent genename %r\n" % parent_gn
                              + "  - genename        %r\n" % genename
                              + "  - event           %r" % event)
                        #raise

                # Nothing gets inserted if it is a duplication.
                insert_missing_spe(parent_sp, ancestor, genename, parent_gn,
                                   node, child, diclinks, ages, event=event)


def search_by_ancestorlist(tree, parse_species_genename, descendants,
                           latest_ancestor=False):
    """Return an iterator over all basal nodes belonging to an ancestor in 
    `descendants`. If `latest_ancestor` is True, return the *speciation* nodes
    where this clade begins.
    """
    def is_in_descendants(node):
        taxon, _ = parse_species_genename(node)
        return taxon if taxon in descendants else None

    if not latest_ancestor:
        #print("")
        stop_at_clade = is_in_descendants
    else:
        def stop_at_clade(node):
            matched_anc = is_in_descendants(node)
            if matched_anc is not None:
                # test for any
                return node.is_leaf() or \
                       not all(is_in_descendants(ch)==matched_anc for ch in node.children)
            else:
                return False
    
    return tree.iter_leaves(is_leaf_fn=stop_at_clade)


def with_dup(leafnames):
    leafspecies = [ultimate_seq2sp(leaf) for leaf in leafnames]
    return (len(leafspecies) > len(set(leafspecies)))


#def get_mostdivergent(nodes, maxsize)
def get_basal(nodes, maxsize):  # ~~> dendro.
    """Identify `maxsize` most basal nodes from a list of sister nodes.

    Return 2 lists:
        - selected basal nodes,
        - excluded nodes (to detach)
    
    "Basal" means: divergence closer to the root.
    
    - `nodes`: a *list* of TreeNode instances;
    - `maxsize`: integer.

    Return ([],[]) if maxsize >= number of leaves.

    sortkey: the tree will be descended following this order.
    An alternative sortkey could be 
        >>> sortkey = lambda n: n.get_closest_leaf(topology_only=False)[1]
    """
    
    #get_dist = lambda n: n.dist
    nodedists = [(node, node.dist) for node in sorted(nodes,key=lambda n:n.dist)]

    # Not enough leaves
    if sum(len(n) for n,_ in nodedists) <= maxsize:  #minsize
        #return [], []
        return nodedists, []
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

    # Sort by distance from the original root:
    sortkey = lambda nodedist: nodedist[1]

    kept = []
    while len(kept) + len(nodedists) < maxsize:
        try:
            node, dist = nodedists.pop(0)  # Descend into the closest divergence.
        except IndexError:
            break

        nextnodes = node.children
        if nextnodes:
            nodedists.extend((nn, nn.dist + dist) for nn in nextnodes)
            nodedists.sort(key=sortkey)
        else:
            kept.append(node)

    #return (to keep, to detach)
    descended_nodes = kept + [n for n,_ in nodedists]
    return descended_nodes[:maxsize], descended_nodes[maxsize:]


def get_data_ete3(tree, nodedist):
    return [(child, child.dist) for child in nodedist[0].children]

#def is_leaf(node):
#    return node.is_leaf()

def reroot_with_outgroup(node, maxsize=0, minsize=-1,
                         is_allowed_outgroup=None,
                         uniq_allowed_value=False,
                         already_taken=None):  # ~~> dendro.reconciled
    """Goes up the tree (towards the root) until it finds outgroup taxa.
    Return (new formatted root (thinned), number of outgroup leaves).
    
    - Only keep at most `maxsize` leaves in the outgroup.
    - Keep all leaves if maxsize < 0.
    - Return None if there is no outgroup.

    If the common ancestor with the outgroup is a multifurcation,
    process all of the sister taxa as a single outgroup.

    Select outgroup leaves to keep those that are the most distant to each other
    (having the most basal latest common ancestor).

    Remark: it doesn't go further than the next immediate outgroup. If you
    want to add the two 2 closest outgroups, give a minsize >= 0.

    param: uniq_allowed_value:
    if is_allowed_outgroup returns a not False or None value, only one leaf
    for each value must be returned.
    Example 1: Use it to exclude paralog sequences.
    Example 2: In combination with minsize>=0, select *mandatory* species.
    """
    logger.debug('Rerooting %r:%s...', node,
                 node.get_ascii(show_internal=False, compact=True))
    logger.debug('maxsize=%d, minsize=%d, is_allowed_outgroup=%s, '
                 'uniq_allowed_value=%s, already_taken=%s',
                 maxsize, minsize, is_allowed_outgroup, uniq_allowed_value,
                 already_taken)
    if maxsize == 0:
        # No need of an outgroup.
        return node, 0
    
    root = node
    # Ignore single child nodes.
    while len(root) == len(node):
        if root.is_root():
            logger.debug('Reached the root of the tree!')
            return None, 0
        # Else go to parent node
        node = root
        root = node.up

    orig_outgroups = node.get_sisters()
    realsize = sum(len(sister) for sister in orig_outgroups)

    if maxsize < 0 and minsize<=0:
        # Attention: ignores the value of is_allowed_outgroups.
        #
        # We include the full outgroup subtrees.
        # Store the recursive child indices to find the same node in the copy.
        ingroup_indices = [root.children.index(node)]
        while root.up and maxsize < -1:
            next_outgroups = root.get_sisters()
            if next_outgroups:  # Else, let's skip the single-child node.
                maxsize += 1
                orig_outgroups += next_outgroups
                realsize += sum(len(sister) for sister in next_outgroups)
                ingroup_indices.insert(0, root.up.children.index(root))
            else:
                ingroup_indices.insert(0, 0)
            root = root.up
        root = root.copy()  # Because we can add the `is_outgroup` mark only on a copy.
        ingroup = root  # re-traverse leafwards to annotate the outgroups.
        for i in ingroup_indices:
            for child in ingroup.children[:i] + ingroup.children[i+1:]:
                child.add_feature('is_outgroup', 1)
            ingroup = ingroup.children[i]

        return root, realsize

    if maxsize > 0 or is_allowed_outgroup is not None:
        # Because root will be copied, we set a marker to remember those nodes.
        for outgroup in orig_outgroups:
            outgroup.add_feature('is_outgroup', 1)
        
        # MAKE A COPY (because this subtree could be reused in another function)
        # NOTE: this does *not* preserve the `up` attribute (parent information).
        orig_root = root
        root = root.copy()

        if uniq_allowed_value:
            assert is_allowed_outgroup is not None, \
                    "`is_allowed_outgroup` is necessary with `uniq_allowed_value`."

            if already_taken is None: already_taken = set((None, False))
            def is_allowed_leaf(node):
                value = is_allowed_outgroup(node)
                if value not in already_taken:
                    already_taken.add(value)
                    return True
                else:
                    return False
        else:
            is_allowed_leaf = (lambda l: True) if is_allowed_outgroup is None else is_allowed_outgroup

        base = [node for node in root.children
                if getattr(node, 'is_outgroup', False)]

        # Now that we know the basal tree:
        # keep only one leaf per identified basal node.

        # Select by proximity
        outgroup_leaves = [(l, d) for basal in base
                           for l,d in iter_distleaves(basal, get_data_ete3)]
        outgroup_leaves.sort(key=lambda leafdist: leafdist[1])
        # Remove unwanted species
        outgroup_leaves = [(l, d) for l,d in outgroup_leaves if is_allowed_leaf(l)][:maxsize]

        realsize = 0
        for basal in base:
            thinned_child = thin(basal.detach(), [l for l,_ in outgroup_leaves])
            if thinned_child is None:
                logger.info('Discard basal node %s without allowed species.',
                            basal.name)
            else:
                #TODO: check the leaf distance. assert outgroup_leaves
                thinned_child.add_feature('is_outgroup', 1)
                root.add_child(child=thinned_child)
                realsize += len(thinned_child)

    for o_outgroup in orig_outgroups:
        o_outgroup.del_feature('is_outgroup')

    while realsize < minsize and root is not None:
        # This block *should* be *inside* this function, because the function
        # itself updates the content of `already_taken`
        logger.debug('Recursive reroot because realsize %d < minsize %d.',
                     realsize, minsize)

        if orig_root.up is None:
            root = None
            break

        if not orig_root.get_sisters():
            # Skip single child nodes.
            orig_root = orig_root.up
            continue
        
        # Get the orig_root.up, copy (with ancestors), replace the orig_root by root.
        orig_root.add_feature('has_ingroup', True)
        nextroot = orig_root.up.copy()  # This doesn't preserve parents.
        orig_root.del_feature('has_ingroup')

        orig_rname = orig_root.get_tree_root().name  # Checking
        #nextrname = nextroot.get_tree_root().name
        #assert nextrname == orig_rname, "Root changed! %s -> %s" % (orig_rname, nextrname)

        nextnode = [node for node in nextroot.children
                    if getattr(node, 'has_ingroup', False)][0]
        orig_root = orig_root.up

        # `root` is the next node.
        # Replace by the thinned version
        if maxsize > 0 or is_allowed_outgroup is not None or uniq_allowed_value:
            nextroot.remove_child(nextnode)
            nextroot.add_child(root)
            nextnode = root

        root, nextsize = reroot_with_outgroup(nextnode,
                                              maxsize-realsize,
                                              0,
                                              is_allowed_outgroup,
                                              uniq_allowed_value,
                                              already_taken)
        realsize += nextsize

    return root, realsize


def save_subtrees(treenb, treefile, ancestor_descendants, ancestor_regexes, #ancgene2sp,
        parse_species_genename,
        diclinks, #treebest=False,
        ages=None, fix_suffix=True, force_mrca=False,
        latest_ancestor=False, #ensembl_version=ENSEMBL_VERSION,
        outdir='.',
        only_dup=False, one_leaf=False, outgroups=0, allowed_outgroups=None,
        reverse=False,
        dry_run=False):
    #print_if_verbose("* treefile: " + treefile)
    #print("treebest = %s" % treebest, file=stderr)
    if isinstance(treefile, ete3.TreeNode):
        tree = treefile
        treefilename = tree.name
    else:
        treefilename = treefile
    
    logger.debug('save_subtrees(\n'
        '   treenb = %d,\n' % treenb +
        '   treefilename = %r,\n' % (treefilename if len(treefilename)<100
                        else (treefilename[:45] + ' ... ' + treefilename[-45:])) +
        '   ancestor_descendants = %.60r [...],\n' % ancestor_descendants  +
        '   ancestor_regexes = %.60r [...],\n' % ancestor_regexes +
        '   parse_species_genename = %r,\n' % parse_species_genename +
        '   diclinks = %.60s [...],\n' % diclinks +
        '   ages = %.60s [...],\n' % ages +
        '   fix_suffix = %r,\n' % fix_suffix +
        '   force_mrca = %r,\n' % force_mrca +
        '   latest_ancestor = %r,\n' % latest_ancestor +
        '   outdir = %r,\n' % outdir +
        '   only_dup = %r,\n' % only_dup +
        '   one_leaf = %r,\n' % one_leaf +
        '   outgroups = %r,\n' % outgroups +
        '   allowed_outgroups = %.50r [...],\n' % allowed_outgroups +
        '   reverse = %r,\n' % reverse +
        '   dry_run = %r)' % dry_run)
    outtrees_set = set() # check whether I write twice the same tree
    outnodes_set = set()

    if not isinstance(treefile, ete3.TreeNode):
        #FIXME: should not allow stdin here:
        if treefile == '-': treefile = '/dev/stdin'
        #try:
        tree, *extratrees = iter_from_ete3(treefile, format=1)
        #except ete3.parser.newick.NewickError as err:
        #    err.args = (err.args[0] + 'ERROR with treefile %r ...' % treefile[:50],)
        #    raise
    else:
        treefile = treefilename
        extratrees = []

    if allowed_outgroups:
        def is_allowed_outgroup(leaf):
            taxon_name = parse_species_genename(leaf)[0]
            if leaf.is_leaf() and taxon_name in allowed_outgroups:
                return taxon_name
            else:
                return False
    else:
        is_allowed_outgroup = None  # This is a valid `is_leaf_fn` argument value in Ete3.

    insert_species_nodes_back(tree, parse_species_genename, diclinks, ages,
                              fix_suffix, force_mrca)
    
    # Output all current features.
    output_features = set.union(set(('is_outgroup',)),
                                *(set(n.features) for n in tree.traverse()))\
                      .difference(('name', 'dist', 'support', 'event'))
    print_if_verbose("* Searching for ancestors:")
    logger.debug('Params: latest_ancestor=%s; only_dup=%s; one_leaf=%s; outgroups=%s',
                 latest_ancestor, only_dup, one_leaf, outgroups)
    for ancestor, descendants in ancestor_descendants.items():
        print_if_verbose(ancestor)
        root_regex, root_repl = ancestor_regexes[ancestor]
        for ancestornodeid, node in enumerate(search_by_ancestorlist(tree,
                                                        parse_species_genename,
                                                        descendants,
                                                        latest_ancestor)):
        #for ancestornodeid, node, root in enumerate(search_by_ancestorlist(tree, ancestorlist, 
        #                                                    latest_ancestor)):
            leafnames = node.get_leaf_names()
            #print(node.name)
            #print(node.get_ascii())
            logger.debug('Got %s', node.name)
            if len(leafnames) > 1 or one_leaf or reverse:
                # check that there is at least one duplication
                if not only_dup or with_dup(leafnames) or reverse:  #FIXME: DaFuQ, this misses dubious dup.
                    ### TODO: when you *know* there can be duplicated node
                    ###       names, use a custom unique id for each node.
                    outname = root_regex.sub(root_repl, node.name).format(
                                    ancestor=ancestor,
                                    rootname=node.name,
                                    nodeid=ancestornodeid,
                                    basename=op.basename(treefile),
                                    treenb=treenb)
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
                    if outname in outtrees_set:
                        #raise RuntimeError("Cannot output twice the same tree (%s)"
                        #                   % outname)
                        logger.error("Should not output twice the same subtree (%s)" % outname)
                    if reverse:
                        # Collapse the node
                        #node.children = []  # Might be better to child.detach()
                        #node.name = outname
                        outnodes_set.add((node, outname))
                        continue
                    outfile = None if outdir == '-' else op.join(outdir, outname + '.nwk')
                    #elif outfile in outfiles_set:
                        # Not sure this case can happen, but better prevent it.
                        # Months later: this case happened.
                        #logger.error("Cannot output twice to %r", outfile)
                    #    raise FileExistsError("Cannot output twice to %r" % outfile)
                    
                    # Add requested outgroups:
                    root, outsize = reroot_with_outgroup(node, outgroups, outgroups,
                                                         is_allowed_outgroup)
                    if not root:
                        logger.warning("No outgroup available for node %r",
                                       node.name)
                        # NOT OUTPUTTED!
                        continue

                    outtrees_set.add(outname)
                    if dry_run:
                        print("node %r (%s) -> %s"
                              % (node.name, ancestor, outfile))
                        root.write(format=1,
                                   format_root_node=True,
                                   outfile=None,
                                   features=output_features)
                    else:
                        #print_if_verbose("Writing to %r." % outfile)
                        logger.info("Writing to %r.", outfile)
                        outtree = root.write(format=1,
                                             format_root_node=True,
                                             outfile=outfile,
                                             features=output_features)
                        if outtree is not None:
                            print(outtree)
                    #outfiles_set.add(outfile)
            else:
                logger.debug('Not outputted')
    #print_if_verbose("Output %d trees." % len(outfiles_set))
    if reverse:
        # Better to collapse nodes after search, because of the parse_species_genename warning.
        for outnode, outname in outnodes_set:
            outnode.children = []
            outnode.name = outname
        root_regex, root_repl = ancestor_regexes['__root__']
        outname = root_regex.sub(root_repl, tree.name).format(
                        rootname=tree.name,
                        nodeid=-1,
                        basename=op.basename(treefile),
                        treenb=treenb)
        outfile = None if outdir == '-' else op.join(outdir, outname + '.nwk')
        if dry_run:
            print("node %r (root) -> %s"
                  % (tree.name, outfile))
            tree.write(format=1,
                       format_root_node=True,
                       outfile=None,
                       features=output_features)
        else:
            #print_if_verbose("Writing to %r." % outfile)
            logger.info("Writing to %r.", outfile)
            outtree = tree.write(format=1,
                                 format_root_node=True,
                                 outfile=outfile,
                                 features=output_features)
            if outtree is not None:
                print(outtree)
        #outfiles_set.add(outfile)
        outtrees_set = set((outname,))

    if extratrees:
        logger.warning('Processing additional trees from: %s', op.basename(treefile))
        for extratree in extratrees:
            extra_outtrees = save_subtrees(treenb, extratree, ancestor_descendants, ancestor_regexes, #ancgene2sp,
                parse_species_genename,
                diclinks, #treebest=False,
                ages=ages, fix_suffix=fix_suffix, force_mrca=force_mrca,
                latest_ancestor=latest_ancestor, #ensembl_version=ENSEMBL_VERSION,
                outdir=outdir,
                only_dup=only_dup, one_leaf=one_leaf, outgroups=outgroups, allowed_outgroups=allowed_outgroups,
                reverse=reverse,
                dry_run=dry_run)
            assert not extra_outtrees & outtrees_set, "Duplicated outtrees from a multi-newick: %s" % (extra_outtrees & outtrees_set)
            outtrees_set |= extra_outtrees
    return outtrees_set


def save_subtrees_process(params, catch_stdout=True):
    logger.info("* Input tree %d: '%s%s'", params[0],
                '...' if len(params[1])>80 else '',
                params[1][-60:])
    ignore_errors = params.pop()
    # Setup the stdout replacement for this subprocess
    #global stdout
    if catch_stdout:
        process_stdout = StringIO()
        #params += (process_stdout,)
        sys.stdout = process_stdout
        get_stdout = process_stdout.getvalue
    else:
        get_stdout = lambda: ''
    try:
        outtrees = save_subtrees(*params)
    except BaseException as err:
        if ignore_errors and not isinstance(err, KeyboardInterrupt):
            logger.info("Ignore %d: %r: %r", params[0], params[1], err)
            return 0, set(), get_stdout()
        else:
            # Allow the loop to finish and compute the summary.
            logger.exception(str(err))
            raise StopIteration()  ## NOT good as the main loop will see this as a normal return...

    # return a value to let Pool.map count the number of results.
    if catch_stdout: sys.stdout = sys.__stdout__
    return 1, outtrees, get_stdout()


def parallel_save_subtrees(treefiles, ancestors, ncores=1, outdir='.',
                           outsub=None, treebest=False, only_dup=False,
                           one_leaf=False, fix_suffix=True, force_mrca=False,
                           latest_ancestor=False, outgroups=0, dry_run=False,
                           ignore_errors=False,
                           ensembl_version=ENSEMBL_VERSION,
                           rename_root_subst=RENAME_ROOT_SUBST, reverse=False):
    ### WARNING: uses global variables here, that are changed by command line
    phyltree = PhylTree.PhylogeneticTree(PHYLTREE_FMT.format(ensembl_version))
    # Crucial point: the pattern alternatives must be sorted by age, so that
    # you don't match an ancestor whose name is contained in its descendant name
    # (like theria is contained in eutheria)
    ancestor_descendants = {}  # Lists of descendants of each given ancestor.
    ancestor_regexes = {}
    rename_root_pattern, rename_root_repl = parse_substitution_expr(rename_root_subst)
    for anc_lowercase in ancestors:
        ancestor = anc_lowercase.capitalize()
        descendants = sorted(phyltree.allDescendants[ancestor],
                             key=lambda anc: -phyltree.ages[anc])
        # Put oldest descendants first.
        ancestor_descendants[anc_lowercase] = descendants
        #FIXME: should be customized when --treebest is set.
        ancestor_regexes[anc_lowercase] = (
            re.compile(rename_root_pattern.format(
                       descendants=r'|'.join(descendants).replace(' ', r'.'),
                       ancestor=ancestor)),
            rename_root_repl)
    if reverse:
        descendants = sorted(phyltree.allNames, key=lambda anc: -phyltree.ages[anc])
        # By default, substitute any taxon name with the PhylTree root.
        ancestor_regexes['__root__'] = (
            re.compile(rename_root_pattern.format(
                       descendants=r'|'.join(phyltree.allNames).replace(' ', r'.'),
                       ancestor=r'|'.join(ancestors)).replace(' ', r'.')),
            rename_root_repl.replace(r'{ancestor}', phyltree.root.replace(' ', r'.')))

    try:
        diclinks = phyltree.dicLinks.common_names_mapper_2_dict()
        ages = phyltree.ages.common_names_mapper_2_dict()
    except AttributeError:
        # Changed myPhylTree so that those objects can be pickled.
        diclinks = phyltree.dicLinks
        ages = phyltree.ages

    # Because of multiprocessing
    global get_species, split_ancestor

    # ~~> dendro.reconciled
    if treebest:
        print_if_verbose("  Reading from TreeBest reconciliation format (S tag and leaf labels '{gene}_{species}')")
        #FIXME: for other applications (TreeRecs, GeneRax, ALE) it's '{species}_{gene}'

        def get_species(node):
            return node.S.replace('.', ' '), node.name.split('_')[0]
        def split_ancestor(node):
            return node.S.replace('.', ' '), node.name
    else:
        ancgene2sp = re.compile(r'('
                            + r'|'.join(list(phyltree.listSpecies) + 
                                        sorted(phyltree.listAncestr,
                                               key=lambda a:len(a),
                                               reverse=True)).replace(' ','\.')
                            + r')([^a-z].*|)$')

        def get_species(node):
            return ultimate_seq2sp(node.name, ensembl_version), node.name
        def split_ancestor(node):
            return split_species_gene(node.name, ancgene2sp)

    this_parse_species_genename = partial(parse_species_genename,
                                          get_species=get_species,
                                          split_ancestor=split_ancestor)

    try:
        outgroups = int(outgroups)
        allowed_outgroups = phyltree.lstEspFull  # Is this really needed?
    except ValueError:
        allowed_outgroups = outgroups.split(',')
        outgroups = len(allowed_outgroups)

    if outsub:
        def format_outdir(treefile):
            return outdir.format(op.basename(treefile).rsplit(outsub, 1)[0])

    else:
        def format_outdir(treefile):
            return outdir.format(op.splitext(op.basename(treefile))[0])

    # NOTE: each arg should be a *list* (because need the .pop() method),
    #       and `ignore_errors` should be the last arg.
    generate_args = [[i, treefile,
                      ancestor_descendants,
                      ancestor_regexes,
                      this_parse_species_genename,
                      #ancgene2sp,
                      diclinks,
                      #treebest,
                      ages,
                      fix_suffix,
                      force_mrca,
                      latest_ancestor,
                      #ensembl_version,
                      format_outdir(treefile),
                      only_dup,
                      one_leaf,
                      outgroups,
                      allowed_outgroups,
                      reverse,
                      dry_run,
                      ignore_errors] for i, treefile in enumerate(treefiles)]

    n_input = len(treefiles)
    logger.info("To process: %d input trees", n_input)
    if ncores > 1:
        install_mp_handler(logger)
        def iter_outputs():
            chunksize = min(20, n_input//ncores + 1)
            with mp.Pool(ncores) as pool:
                yield from pool.imap_unordered(save_subtrees_process,
                                               generate_args,
                                               chunksize)
    else:
        def iter_outputs():
            return (save_subtrees_process(args, catch_stdout=False) for args in generate_args)

    progress = 0
    all_outtrees = set()
    matching_inputs = 0
    for return_value, outtrees, stdout_txt in iter_outputs():
        print(stdout_txt, end='')
        progress += return_value
        dup_outtrees = all_outtrees & outtrees
        assert not dup_outtrees, \
            "Duplicated outtrees from different input trees: %s" % dup_outtrees
        all_outtrees |= outtrees
        matching_inputs += bool(outtrees)  # +1 if not empty.

    logger.info("  \nFinished processing %d/%d input trees.\n"
                "  Output %d trees found in %d input trees%s%s.",
                progress, n_input, len(all_outtrees), matching_inputs,
                ' [DRY-RUN]' if dry_run else '',
                ' [ignore errors]' if ignore_errors else '')
    return all_outtrees


def parse_treelist(treelist_file):
    try:
        if treelist_file=='-':
            return [line.rstrip() for line in stdin]
        with open(treefiles_file) as stream:
            return [line.rstrip() for line in stream]
    except BrokenPipeError:
        logger.warning('Broken pipe.')


if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)s:%(name)s:%(module)s:%(funcName)s:%(message)s')
    logger.setLevel(logging.INFO)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("treefile")
    parser.add_argument("ancestors", nargs='+')
    intype_parser = parser.add_mutually_exclusive_group()
    intype_parser.add_argument("-f", "--fromfile", action="store_true",
                               help="read treefile"\
                                    " names from the first argument")
    intype_parser.add_argument("-m", "--multi-newick", action="store_true",
                               help="Parallelize over multiple trees from a single input file.")
    parser.add_argument("-o", "--outdir", default='./{0}',
                        help="'-' to output trees to stdout. [%(default)s]. "\
                             "The basename is automatic based on the root node"\
                             " (but you can use --rename-root-subst).")
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
    parser.add_argument("-O", "--outgroups", metavar='N', default=0, 
            help="Save the subtree including outgroups: `%(metavar)s` can be: "\
                "1) a positive integer: the exact number of outgroup *leaves*; "\
                "2) a negative number: the number of *ancestors* to go back "\
                "(return complete sister clades); "\
                             "3) a comma-separated list of species to "\
                             "restrict to a given set of species."
                             "NOTE: with `-l`, the outgroup can be a paralog "\
                             "tree, otherwise it's necessarily in sister "\
                             "species. If specified, DO NOT output trees "\
                             "without outgroups. [Default $(metavar)s=%(default)s]")
    parser.add_argument("-r", "--rename-root-subst", "--rr",
                        default=RENAME_ROOT_SUBST,
                        help="sed substitution expression to build output "\
                            "filenames from the root node information." \
                            "ancestor_regexes. Customise in order to avoid " \
                            "duplicated outnames. [%(default)r].\n"\
                            "Available fields are: {descendants} {ancestor} "\
                            "(replaced) {rootname} {basename} {treenb} {nodeid}"\
                            "(replacement).\n" \
                            "Example to always prefix with the input tree index: "\
                            "'s/^/{treenb}/'.")
    parser.add_argument("-R", "--reverse", action="store_true",
                        help="Output the super tree (ancestral) with ancestors as leaves.")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="only print out the output files it would produce")
    parser.add_argument("--ncores", type=int, default=1, help="Number of cores")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument("-i", "--ignore-errors", action="store_true", 
                        help="On error, print the error and continue the loop.")
    #parser.add_argument("-O", "--to-stdout", action="store_true")
    
    args = parser.parse_args()
    dargs = vars(args)

    PHYLTREE_FMT = dargs.pop("phyltree_fmt")
    #ENSEMBL_VERSION = dargs.pop("ensembl_version")
    #ANCGENE_START = dargs.pop("ancgene_start")
    #ANCGENE2SP = re.compile(ANCGENE2SP_PATTERN % ANCGENE_START)

    verbosity = dargs.pop("verbose")
    if not verbosity:
        # mute the custom print function
        def print_if_verbose(*args, **kwargs):
            pass
    elif verbosity > 1:
        logger.setLevel(logging.DEBUG)

    fromfile = dargs.pop('fromfile')
    multi_newick = dargs.pop("multi_newick")
    if fromfile:
        treefiles = parse_treelist(dargs.pop("treefile"))
    elif multi_newick:
        treefile = dargs.pop("treefile")
        try:
            with (open(treefile) if treefile != '-' else stdin) as tree_input:
            #    treefiles = [txt + ';' for txt in tree_input.read().split(';')]
                treefiles = list(read_multinewick(tree_input))
        except BrokenPipeError:
            logger.warning('Broken pipe.')
        dargs['outdir'] = dargs['outdir'].format(op.basename(treefile))
    else:
        treefiles = [dargs.pop("treefile")]

    setrecursionlimit(50000)
    parallel_save_subtrees(treefiles, **dargs)

