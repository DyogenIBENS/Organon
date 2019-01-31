#!/usr/bin/env python3
# -*- coding:utf8 -*-

"""Interactively display a tree using ete3."""

from sys import stdin, stderr
import os.path as op
import ete3
import argparse
#import fileinput


def fileinputtrees():
    newicks = ['']
    #lines = [stdin.readline()]

    #while lines:
    for line in stdin:

        line_chunks = line.rstrip().split(';')
        #[chunk for chunk in line.rstrip().split(';') if chunk]
        
        newicks[-1] += line_chunks.pop(0)
        for chunk in line_chunks:
            newicks[-1] += ';'
            newicks.append(chunk)  # Works even when chunk == ''

    if not newicks[-1]:
        newicks.pop()
    return newicks

    ## Méthode yolo
    #return [newick.rstrip() + ';' for newick in stdin.read().split(';')
    #        if newick.rstrip()]

def make_treename(newickfile):
    try:
        newickfile = newickfile.name
    except AttributeError:
        pass

    if op.isfile(newickfile):
        return op.splitext(op.basename(newickfile))[0]


def main(newickfiles, subnewick_format=1, rootnode=None, ladderize=False,
         show_internal=False, show_dists=False, nodesize=1, text=False,
         compact=False, output=None, parser='ete3'):
    #newicks = fileinput.input(files=newickfiles)
    if parser.lower() == 'ete3':
        def parse_trees(newickfiles):
            if not newickfiles:
                newickfiles = fileinputtrees()
                #print(newickfiles)
            for newick in newickfiles:
                #print('subnewick format = %d' % subnewick_format)
                try:
                    tree = ete3.Tree(newick, format=subnewick_format)
                except:
                    print("newick=%r" % newick, file=stderr)
                    raise
                yield tree, make_treename(newick)
    else:
        from dendron.parsers import parserchoice
        from dendron.converters import converterchoice
        treeparser = parserchoice[parser]
        treeconverter = converterchoice[parser]['ete3']
        def parse_trees(treefiles):
            if not treefiles:
                treefiles = [stdin]
            for treefile in treefiles:
                for tree in treeparser(treefile):
                    yield treeconverter(tree), make_treename(treefile)

    for i, (tree, treename) in enumerate(parse_trees(newickfiles), start=1):
        print('Tree %d' % i)
        display_onetree(tree, treename, rootnode, ladderize, show_internal,
                        show_dists, nodesize, text, compact, output)


def display_onetree(tree, treename=None, rootnode=None, ladderize=False,
                    show_internal=False, show_dists=False, nodesize=1,
                    text=False, compact=False, output=None):

    if rootnode:
        #print(tree.search_nodes(name=rootnode))
        tree = tree.search_nodes(name=rootnode)[0]

    if ladderize:
        tree.ladderize()

    if text:
        show_attr=['dist'] if show_dists else None
        print(tree.get_ascii(show_internal=show_internal, compact=compact,
                             attibutes=show_attr))
    else:
        mynodestyle = ete3.NodeStyle(size=nodesize)

        ts = ete3.TreeStyle()
        if show_dists:
            ts.show_branch_length = True
        
        def mybasiclayout(node):
            node.set_style(mynodestyle)

        if show_internal:
            nameface = ete3.faces.AttrFace("name")
            def mylayout(node):
                mybasiclayout(node)
                if not node.is_leaf():
                    ete3.faces.add_face_to_node(nameface, node, column=0)
        else:
            mylayout = mybasiclayout

        if output:
            tree.render(output, tree_style=ts, layout=mylayout)
        else:
            tree.show(tree_style=ts, layout=mylayout, name=treename)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('newickfiles', nargs='*')
    parser.add_argument('-f', dest='subnewick_format', default=1, type=int,
                        help='subnewick format (from 0 to 9, or 100), ' \
                             'see ete3 documentation. [%(default)s]')
    parser.add_argument('-r', dest='rootnode', help='rootnode: name of an ' \
                        'internal node to take as root.')
    parser.add_argument('-l', dest='ladderize', action='store_true',
                        help='Ladderize tree before display')
    parser.add_argument('-i', dest='show_internal', action='store_true', 
                        help='Display internal nodes names')
    parser.add_argument('-d', dest='show_dists', action='store_true', 
                        help='Display branch lengths')
    parser.add_argument('-s', '--nodesize', type=int, default=1)
    parser.add_argument('-t', '--text', action='store_true', 
                        help='Print as text')
    parser.add_argument('-c', '--compact', action='store_true',
                        help='Compact output (only with --text)')
    parser.add_argument('-o', '--output', help='save to output file')
    parser.add_argument('-p', '--parser', default='ete3',
                        choices=['ete3', 'prottree', 'phyltree'],
                        help='Implementation to use. [%(default)s]')
    
    args = parser.parse_args()
    
    main(**vars(args))

