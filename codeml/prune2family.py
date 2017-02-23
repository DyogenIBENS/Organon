#!/usr/bin/env python3

from __future__ import print_function

import re
import sys
import os.path
import argparse
import multiprocessing as mp

import ete3
import LibsDyogen.myPhylTree as PhylTree
from select_leaves_from_specieslist import convert_gene2species


ENSEMBL_VERSION = 85
PHYLTREE_FMT = "/users/ldog/alouis/ws2/GENOMICUS_SVN/data{0}/PhylTree.Ensembl.{0}.conf"
NEW_DUP_SUFFIX = re.compile(r'\.[A-Za-z`]+$')

#SPLIT_SPECIES_GENE = re.compile()

def print_if_verbose(*args, **kwargs):
    print(*args, **kwargs)


def split_species_gene(nodename):
    """Split a node name into its two parts (taxon + ancestral gene name)."""
    try:
        idx = nodename.index('ENS')
    except ValueError:
        try:
            idx = nodename.index('FBgn') # Drosophila
        except ValueError:
            try:
                idx = nodename.index('WBGene') # Caenorhabditis
            except ValueError:
                try:
                    idx = nodename.index('Y')
                except ValueError:
                    try:
                        idx = nodename.index('Q0')
                    except ValueError as err:
                        err.args += ("ERROR: Invalid nodename %r" % nodename,)
                        raise
    return nodename[:idx].replace('.', ' '), nodename[idx:]


def name_missing_links(parent_sp, ancestor, genename, parent_node_name,
                       child_name, diclinks, ages=None):
    """given two taxa and the child gene name, return each missing node
    inbetween."""
    try:
        ancestor_lineage = diclinks[parent_sp][ancestor]
    except KeyError as err:
        err.args += ("child node : %s (%r)" % (child_name, ancestor),
                     "parent node: %s (%r)" % (parent_node_name, parent_sp))
        raise err
    # Links from diclinks must be corrected: unnecessary nodes removed, i.e
    # nodes with a single child and age 0. (introducing potential errors)
    # BUT: do not remove species (age = 0) !
    if ages:
        links = [link for link in ancestor_lineage[1:-1] if ages[link] > 0]

    # If doesn't match species tree
    # same ancestor as child is possible for duplication node
    # So check for length 1 or 2
    new_node_names = []
    if links:
        new_node_names = [link + genename for link in links]

    if ages:
        new_branch_dists = []
        total_len = ages[ancestor_lineage[-1]] - ages[ancestor_lineage[0]]
        for link_parent, link in zip(ancestor_lineage[:-1],
                                     ancestor_lineage[ 1:]):
            new_dist = float(ages[link] - ages[link_parent])
            try:
                new_branch_dists.append(new_dist / total_len)
            except ZeroDivisionError as err:
                err.args += (parent_sp, ancestor, "please check ages.")
                print("WARNING:", err, file=sys.stderr)
                new_branch_dists.append(0)
    else:
        total_len = len(ancestor_lineage) - 1
        new_branch_dists = [1./total_len] * total_len

    return new_node_names, new_branch_dists
    # TODO: compute intermediate distances proportionally to the age of each
    # intermediate taxon. Will however be incorrect if one node is a duplication.


def insert_nodes(new_node_names, parent, child, new_dist_ratios=None):
    # conserve original distance
    if not new_dist_ratios:
        # split equally
        n_new_branches = len(new_node_names) + 1
        new_branch_dists = [float(child.dist) / n_new_branches] * n_new_branches
    else:
        new_branch_dists = [float(child.dist) * r for r in new_dist_ratios]

    new_node = parent
    print_if_verbose("Inserted nodes: ", end=" ")
    for new_name, new_dist in zip(new_node_names, new_branch_dists[:-1]):
        new_node = new_node.add_child(name=new_name, dist=new_dist)
        new_node.add_feature('reinserted', True)
        print_if_verbose(new_name, end=" ")
    print_if_verbose()
    # Move the child on top of the created intermediate links
    new_node.add_child(child=child.detach(), dist=new_branch_dists[-1])


def insert_missing_links(parent_sp, ancestor, genename, parent_node, child,
                         diclinks, ages=None):
    new_node_names, new_branch_dists = name_missing_links(parent_sp, ancestor,
                                                          genename,
                                                          parent_node.name,
                                                          child.name,
                                                          diclinks,
                                                          ages)
    if new_node_names:
        print_if_verbose("Nodes to insert: %s" % new_node_names)
        insert_nodes(new_node_names, parent_node, child, new_branch_dists)
        #return True
    # return false if no node was added
    #return False


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
            ancestor, genename = split_species_gene(node.name)
        
        parent_node = node.up
        parent_ancestor, _ = split_species_gene(parent_node.name)

        insert_missing_links(parent_ancestor, ancestor, genename, parent_node,
                             node, diclinks, ages)


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


def insert_species_nodes_back(tree, diclinks, ages=None):
    print_if_verbose("* Insert missing nodes:")
    ### Insert childs *while* iterating.
    ### Beware not to iterate on these new nodes:
    ### iterate from leaves to root ('postorder') and insert between node and child
    for node in tree.traverse('postorder'):
        print_if_verbose(" " + node.name)
        if node.children:
            parent_sp, parent_gn = split_species_gene(node.name)
            child_sp = []
            child_gn = []
            for child in node.children:
                print_if_verbose("  - child %r" % child.name)
                if child.is_leaf():
                    ancestor = convert_gene2species(child.name)
                    genename = child.name
                else:
                    ancestor, genename = split_species_gene(child.name)
                child_sp.append(ancestor)
                child_gn.append(genename)

            test_same_species = [parent_sp == sp for sp in child_sp]
            if all(test_same_species):
                # duplication.
                # check suffixes
                event = 'dup'
            elif not any(test_same_species):
                # speciation
                # check suffixes
                event = 'spe'
            else:
                print_if_verbose("implicit dup+spe")
                # there is a dup + a spe. Need to add nodes
                added_suffixes = [suffix_list(parent_gn, gn) for gn in child_gn]
                added_suffixes = [suf[0].lstrip('.') if suf else None for suf in added_suffixes]
                suff_count = 1 #len(set((suf for suf in added_suffixes if suf)))
                for i, (child, ancestor, genename, suf) in \
                        enumerate(zip(node.children, child_sp, child_gn,
                                        added_suffixes)):
                    if ancestor != parent_sp:
                        # append proper suffix if there isn't one.
                        if not suf:
                            suf = chr(96+suff_count)
                            while suf in added_suffixes:
                                suff_count += 1
                                suf = chr(96+suff_count)
                            added_suffixes.append(suf)
                            genename = parent_gn + '.' + suf
                            child_gn[i] = genename
                        # Insert back the needed speciation event
                        spe_node_name = parent_sp + genename
                        print_if_verbose("Inserted node (spe-dup):",
                                         spe_node_name)
                        # Then check the succession of ancestors
                        new_node_names, new_branch_dists = name_missing_links(
                                                                parent_sp,
                                                                ancestor,
                                                                genename,
                                                                spe_node_name,
                                                                child.name,
                                                                diclinks,
                                                                ages)
                        new_node_names.insert(0, spe_node_name)
                        # adjust the length of this dup-spe branch
                        n_new_nodes = len(new_node_names)
                        new_branch_dists = [1./(n_new_nodes+1)] + \
                                            [r * n_new_nodes/(n_new_nodes+1) \
                                                for r in new_branch_dists]
                        insert_nodes(new_node_names, node, child, new_branch_dists)
                        child_sp[i] = parent_sp
                event = 'dup'
            
            for child, ancestor, genename in zip(node.children, child_sp, child_gn):
                if not child.is_leaf():
                    try:
                        assert suffixes_ok(parent_gn, genename, event)
                    except AssertionError:
                        print("WARNING: inconsistent suffixes:\n",
                              "  - node            %r\n" % node.name,
                              "  - children        %s\n" % (node.children,),
                              "  - parent sp       %r\n" % parent_sp,
                              "  - species         %r\n" % ancestor,
                              "  - parent genename %r\n" % parent_gn,
                              "  - genename        %r\n" % genename,
                              "  - event           %r" % event, file=sys.stderr)
                        #raise

                insert_missing_links(parent_sp, ancestor, genename, node,
                                     child, diclinks, ages)


def search_by_ancestorlist(tree, ancestorlist):
    def stop_at_any_ancestor(node):
        return any(node.name.startswith(anc) for anc in ancestorlist)
    return tree.iter_leaves(is_leaf_fn=stop_at_any_ancestor)


def save_subtrees_byspecieslist(tree, specieslist, outdir='.'):
    ancestor = PHYLTREE.lastCommonAncestor(specieslist)
    for node in search_by_ancestorspecies(tree, ancestor):
        outfile = os.path.join(outdir, node.name + '.nwk')
        node.write(format=1, outfile=outfile)


def with_dup(leafnames):
    leafspecies = [convert_gene2species(leaf) for leaf in leafnames]
    return (len(leafspecies) > len(set(leafspecies)))


def save_subtrees(treefile, ancestorlists, ancestor_regexes, diclinks,
                  ages=None, outdir='.', only_dup=False, dry_run=False):
    #print_if_verbose("* treefile: " + treefile)
    outfiles_set = set() # check whether I write twice to the same outfile
    try:
        tree = ete3.Tree(treefile, format=1)
    except ete3.parser.newick.NewickError as err:
        err.args += ('ERROR with treefile %r' % treefile,)
        raise
    insert_species_nodes_back(tree, diclinks, ages)
    print_if_verbose("* Searching for ancestors:")
    for ancestor, ancestorlist in ancestorlists.items():
        print_if_verbose(ancestor)
        ancestor_regex = ancestor_regexes[ancestor]
        for node in search_by_ancestorlist(tree, ancestorlist):
            leafnames = node.get_leaf_names()
            #print(node.name)
            #print(node.get_ascii())
            if len(leafnames) > 1:
                # check that there is at least one duplication
                if not only_dup or with_dup(leafnames):
                    outname = ancestor_regex.sub(node.name, ancestor)
                    outfile = os.path.join(outdir, node.name + '.nwk')
                    if outfile in outfiles_set:
                        # Not sure this case can happen, but better prevent it
                        print("ERROR: Cannot output twice to %r" % outfile,
                              file=sys.stderr)
                        sys.exit(1)
                    else:
                        outfiles_set.add(outfile)

                    if dry_run:
                        print(outfile)
                    else:
                        print_if_verbose("Writing to %r." % outfile)
                        node.write(format=1,
                                   format_root_node=True,
                                   outfile=outfile,
                                   features=["reinserted"])

def save_subtrees_process(params):
    print("* Input tree: %r" % params[0])
    ignore_errors = params.pop()
    try:
        save_subtrees(*params)
    except BaseException as err:
        if ignore_errors:
            print("Ignore %r: %r" % (params[0], err), file=sys.stderr)
            return 0
        else:
            raise
    return 1 # return a value to let Pool.map count the number of results.


def parallel_save_subtrees(treefiles, ancestors, ncores=1, outdir='.',
                           only_dup=False, dry_run=False, ignore_errors=False):
    phyltree = PhylTree.PhylogeneticTree(PHYLTREE_FMT.format(ENSEMBL_VERSION))
    ancestorlists = {}
    ancestor_regexes = {}
    for anc_lowercase in ancestors:
        ancestor = anc_lowercase.capitalize()
        ancestorlist = sorted(phyltree.allDescendants[ancestor],
                              key=lambda anc: -phyltree.ages[anc])
        ancestorlists[anc_lowercase] = ancestorlist
        ancestor_regexes[anc_lowercase] = re.compile('^(%s)(?=ENSGT)' % \
                                    '|'.join(ancestorlist).replace(' ', '.'))

    diclinks = phyltree.dicLinks.common_names_mapper_2_dict()
    ages = phyltree.ages.common_names_mapper_2_dict()
    generate_args = [[treefile,
                      ancestorlists,
                      ancestor_regexes,
                      diclinks,
                      ages,
                      outdir.format(os.path.splitext(os.path.basename(treefile))[0]),
                      only_dup,
                      dry_run,
                      ignore_errors] for treefile in treefiles]

    n_input = len(treefiles)
    print("To process: %d input trees" % n_input)
    if ncores > 1:
        pool = mp.Pool(ncores)
        return_values = pool.map(save_subtrees_process, generate_args)
        progress = sum(return_values)
    else:
        progress = 0
        for args in generate_args:
            progress += 1
            save_subtrees_process(args)
    print("Finished processing %d/%d input trees" % (progress, n_input))


def parse_treefiles(treefiles_file):
    with open(treefiles_file) as stream:
        return [line.rstrip() for line in stream]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("treefile")
    parser.add_argument("ancestors", nargs='+')
    parser.add_argument("--fromfile", action="store_true", help="read treefile"\
                        " names from the first argument")
    parser.add_argument("-o", "--outdir", default='./{0}', help="[%(default)s]")
    parser.add_argument("--only-dup", action="store_true",
                        help="do not extract trees that don't have at least "\
                             "one duplication")
    parser.add_argument("-e", "--ensembl-version", type=int, default=85,
                        help="[%(default)s]")
    parser.add_argument("-p", "--phyltree-fmt", default=PHYLTREE_FMT,
                        help="Phylogenetic species tree "\
                        "in LibsDyogen PhylTree format. Can contain the string"\
                        " '{0}' which will be replaced by the Ensembl version"\
                        " [%(default)s]")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="only print out the output files it would produce")
    parser.add_argument("--ncores", type=int, default=1, help="Number of cores")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-i", "--ignore-errors", action="store_true", 
                        help="On error, print the error and continue the loop.")
    
    args = parser.parse_args()
    dargs = vars(args)

    PHYLTREE_FMT = dargs.pop("phyltree_fmt")
    ENSEMBL_VERSION = dargs.pop("ensembl_version")

    if not dargs.pop("verbose"):
        # mute the custom print function
        def print_if_verbose(*args, **kwargs):
            pass

    if dargs.pop("fromfile"):
        treefiles = parse_treefiles(dargs.pop("treefile"))
    else:
        treefiles = [dargs.pop("treefile")]
    parallel_save_subtrees(treefiles, **dargs)

