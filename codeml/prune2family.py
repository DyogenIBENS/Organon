#!/usr/bin/env python2.7

import re
import sys
import os.path
import argparse
import multiprocessing as mp

import ete3
import LibsDyogen.utils.myPhylTree as PhylTree
from select_leaves_from_specieslist import convert_gene2species


ENSEMBL_VERSION = 85
PHYLTREE_FMT = "/users/ldog/alouis/ws2/GENOMICUS_SVN/data{0}/PhylTree.Ensembl.{0}.conf"
NEW_DUP_SUFFIX = re.compile(r'\.[A-Za-z`]+$')

def split_species_gene(nodename):
    """When genename is the concatenation of the species and gene names"""
    idx = nodename.index('ENSGT')
    return nodename[:idx].replace('.', ' '), nodename[idx:]


def insert_missing_links(parent_sp, ancestor, genename, parent_node, child,
                         diclinks):
    try:
        ancestor_lineage = diclinks[parent_sp][ancestor] 
    except KeyError:
        print >>sys.stderr, "child node : %s (%r)" % (child.name,
                                                      ancestor)
        print >>sys.stderr, "parent node: %s (%r)" % (node.name,
                                                      parent_sp)
        raise
    # If doesn't match species tree
    # same ancestor as child is possible for duplication node
    # So check for length 1 or 2
    if len(ancestor_lineage) > 2:
        # conserve original distance
        n_new_branches = len(ancestor_lineage) - 1
        dist_new_branches = float(child.dist) / n_new_branches
        # Add missing links
        new_node = parent_node
        for link in ancestor_lineage[1:-1]:
            new_node = new_node.add_child(name=(link + genename),
                                          dist=dist_new_branches)
        # Move the child on top of the created intermediate links
        new_node.add_child(child=child.detach(),
                           dist=dist_new_branches)


def add_species_nodes_back(tree, diclinks):
    """WRONG. DO NOT USE.
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
                             node, diclinks)


def suffixes_ok(parent, child, event):
    """parent: genename
        child: genename
        event: 'dup' or 'spe'"""
    #TODO: what to do with leaves (modern gene names)
    if event == 'spe':
        return parent == child
    elif event =='dup':
        new_suffix = child[(len(parent):]
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
        raise RuntimeError("parent %r and child %r are not in the same lineage")

    difference = child[len(parent):]
    suffixes = []
    match = NEW_DUP_SUFFIX.search(difference)
    while match:
        suffixes.append(match.group())
        difference = NEW_DUP_SUFFIX.sub('', difference)
    return suffixes


def insert_species_nodes_back(tree, diclinks):
    for node in tree.iter_descendants():
        if node.children:
            parent_sp, parent_gn = split_species_gene(node.name)
            child_sp = []
            child_gn = []
            for child in node.children:
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
                # there is a dup + a spe. Need to add nodes
                for child, ancestor, genename in zip(node.children, child_sp,
                                                     child_gn):
                    if ancestor != parent_sp:
                        # Insert back the needed speciation event
                        spe_node = node.add_child(name=(parent_sp + genename))
                        spe_node.add_feature('reinserted', True)
                        spe_node.add_child(child=child.detach())
                        # Then check the succession of ancestors
            
            for child, ancestor, genename in zip(node.children, child_sp, child_gn):
                if not child.is_leaf()
                    assert suffixes_ok(parent_gn, genename, event)
                insert_missing_links(parent_sp, ancestor, genename,
                                     node, child, diclinks)



def search_by_ancestorlist(tree, ancestorlist):
    def stop_at_any_ancestor(node):
        return any(node.name.startswith(anc) for anc in ancestorlist)
    return tree.iter_leaves(is_leaf_fn=stop_at_any_ancestor)


def save_subtrees_byspecieslist(tree, specieslist, outdir='.'):
    ancestor = PHYLTREE.lastCommonAncestor(specieslist)
    for node in search_by_ancestorspecies(tree, ancestor):
        outfile = os.path.join(outdir, node.name + '.nwk')
        node.write(format=1, outfile=outfile)

#def stop_at_duplicates(values):
#    for values in
def with_dup(leafnames):
    leafspecies = [convert_gene2species(leaf) for leaf in leafnames]
    return (len(leafspecies) > len(set(leafspecies)))


def save_subtrees(treefile, ancestorlists, ancestor_regexes, diclinks,
                  outdir='.', only_dup=False, dry_run=False):
    outfiles_set = set() # check whether I write twice to the same outfile
    try:
        tree = ete3.Tree(treefile, format=1)
    except ete3.parser.newick.NewickError as e:
        print >>sys.stderr, 'ERROR with treefile %r' % treefile
        raise
    insert_species_nodes_back(tree, diclinks)
    for ancestor, ancestorlist in ancestorlists.iteritems():
        ancestor_regex = ancestor_regexes[ancestor]
        for node in search_by_ancestorlist(tree, ancestorlist):
            leafnames = node.get_leaf_names()
            #print node.name
            #print node.get_ascii()
            if len(leafnames) > 1:
                # check that there is at least one duplication
                if not only_dup or with_dup(leafnames):
                    outname = ancestor_regex.sub(node.name, ancestor)
                    outfile = os.path.join(outdir, node.name + '.nwk')
                    if outfile in outfiles_set:
                        # Not sure this case can happen, but better prevent it
                        print >>sys.stderr, \
                                "ERROR: Cannot output twice to %r" % outfile
                        sys.exit(1)
                    else:
                        outfiles_set.add(outfile)

                    if dry_run:
                        print outfile
                    else:
                        node.write(format=1,
                                   format_root_node=True,
                                   outfile=outfile)

def save_subtrees_process(params):
    save_subtrees(*params)


def parallel_save_subtrees(treefiles, ancestors, ncores=1, outdir='.',
                           only_dup=False, dry_run=False):
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
    pool = mp.Pool(ncores)

    diclinks = phyltree.dicLinks.common_names_mapper_2_dict()
    generate_args = [(treefile,
                        ancestorlists,
                        ancestor_regexes,
                        diclinks,
                        outdir,
                        only_dup,
                        dry_run) for treefile in treefiles]

    pool.map(save_subtrees_process, generate_args)


def parse_treefiles(treefiles_file):
    with open(treefiles_file) as stream:
        return [line.rstrip() for line in stream]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("treefile")
    parser.add_argument("ancestors", nargs='+')
    parser.add_argument("--fromfile", action="store_true", help="read treefile"\
                        " names from the first argument")
    parser.add_argument("--ncores", type=int, default=1, help="Number of cores")
    parser.add_argument("-o", "--outdir", default='.')
    parser.add_argument("--only-dup", action="store_true",
                        help="do not extract trees that don't have at least "\
                             "one duplication")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="only print out the output files it would produce")
    parser.add_argument("-e", "--ensembl-version", type=int, default=85)
    parser.add_argument("-p", "--phyltree-fmt", default=PHYLTREE_FMT,
                        help="Phylogenetic species tree "\
                        "in LibsDyogen PhylTree format. Can contain the string"\
                        " '{0}' which will be replaced by the Ensembl version")
    
    args = parser.parse_args()
    dargs = vars(args)

    PHYLTREE_FMT = dargs.pop("phyltree_fmt")
    ENSEMBL_VERSION = dargs.pop("ensembl_version")

    if dargs.pop("fromfile"):
        treefiles = parse_treefiles(dargs.pop("treefile"))
    else:
        treefiles = [dargs.pop("treefile")]
    parallel_save_subtrees(treefiles, **dargs)

