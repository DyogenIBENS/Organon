#!/usr/bin/env python3
# -*- coding:utf8 -*-

"""Interactively display a tree (in newick format) using ete3."""

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


def main(newickfiles, subnewick_format=1, rootnode=None, show_internal=False,
         nodesize=1, text=False, compact=False, output=None):
    #newicks = fileinput.input(files=newickfiles)
    if not newickfiles:
        newickfiles = fileinputtrees()
        #print(newickfiles)

    for i, newick in enumerate(newickfiles, start=1):
        print('Tree %d' % i)
        display_onetree(newick, subnewick_format, rootnode, show_internal,
                        nodesize, text, compact, output)


def display_onetree(newick, subnewick_format=1, rootnode=None,
                    show_internal=False, nodesize=1, text=False, compact=False,
                    output=None):
    #print('subnewick format = %d' % subnewick_format)
    try:
        tree = ete3.Tree(newick, format=subnewick_format)
    except:
        print("newick=%r" % newick, file=stderr)
        raise

    if rootnode:
        #print(tree.search_nodes(name=rootnode))
        tree = tree.search_nodes(name=rootnode)[0]

    if text:
        print(tree.get_ascii(show_internal=show_internal, compact=compact))
    else:
        mynodestyle = ete3.NodeStyle(size=nodesize)

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

        treename = op.splitext(op.basename(newick))[0] if \
                        op.isfile(newick) else None
        if output:
            tree.render(output, layout=mylayout)
        else:
            tree.show(layout=mylayout, name=treename)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('newickfiles', nargs='*')
    parser.add_argument('-f', dest='subnewick_format', default=1, type=int,
                        help='subnewick format (from 0 to 9, or 100), ' \
                             'see ete3 documentation. [%(default)s]')
    parser.add_argument('-r', dest='rootnode', help='rootnode: name of an ' \
                        'internal node to take as root.')
    parser.add_argument('-i', dest='show_internal', action='store_true', 
                        help='Display internal nodes names')
    parser.add_argument('-s', '--nodesize', type=int, default=1)
    parser.add_argument('-t', '--text', action='store_true', 
                        help='Print as text')
    parser.add_argument('-c', '--compact', action='store_true',
                        help='Compact output (only with --text)')
    parser.add_argument('-o', '--output', help='save to output file')
    
    args = parser.parse_args()
    
    main(**vars(args))

