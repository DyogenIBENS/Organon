#!/usr/bin/env python3


import sys
import re
from collections import namedtuple
from itertools import zip_longest

from dendron.climber import dfw_descendants_generalized, \
                            dfw_pairs_generalized
from dendron.sorter import ladderize


### TODO: rename as `printer`

## Ascii drawing ##

### branch symbols

# ┏ ┣ ╋ ┫ ┗ ━ ┳ ┻ ┃ = 250f 2523 254b 252b 2517 2501 2533 253b 2503
# ┌ ├ ┼ ┤ └ ─ ┬ ┴ │ = 250c 251c 253c 2524 2514 2500 252c 2534 2502
# top corner, right tack, left tack, bottom corner, branch, empty
CharTuple = namedtuple('chars',
"leaf topleaf bottomleaf node topnode bottomnode hbranch vbranch empty")

chars = CharTuple(
         leaf=       u' \u251c ',
         topleaf=    u' \u250c ',
         bottomleaf= u' \u2514 ',
         node=       u' \u253c ',
         topnode=    u' \u252c ',
         bottomnode= u' \u2534 ',
         hbranch=    u' \u2500 ',
         vbranch=    u' \u2502 ',
         empty=      '   ')

#chars = CharTuple(
#         leaf=       u' \u251c\u2500'     ,   # u' \u2523\u2501'     , 
#         topleaf=    u' \u250c\u2500'     ,   # u' \u250f\u2501'     , 
#         bottomleaf= u' \u2514\u2500'     ,   # u' \u2517\u2501'     , 
#         node=       u'\u2500\u253c\u2500',   # u'\u2501\u254b\u2501', 
#         topnode=    u'\u2500\u252c\u2500',   # u'\u2501\u2533\u2501', 
#         bottomnode= u'\u2500\u2534\u2500',   # u'\u2501\u253b\u2501', 
#         hbranch=    u'\u2500'*3          ,   # u'\u2501'*3          , 
#         vbranch=    u' \u2502 '          ,   # u' \u2503 '          , 
#         empty=      '   ')


def draw_trees(tree, root, get_children=None, fmt_node=None, *fmt_args):
    """Print a tree using ascii characters for branches"""
    positions = {root: (0,0)}
    if not get_children:
        get_children = lambda tree,node: tree.get(node, [])
    if not fmt_node:
        #fmt_node = lambda node, info: "%s (%s)" %(node, info[node].upper()[0])
        fmt_node = lambda node: "%s" % node
    rootlabel = fmt_node(root, *fmt_args)
    drawing = [[rootlabel]] # each inner list is a line
    width = 1
    for g1, g2 in dfw_pairs_generalized(tree, get_children, [(None,root)]):
        g1x, g1y = positions[g1]
        if len(drawing[g1y]) == g1x+1:
            # g1 does not have a descendant on the same line.
            label = fmt_node(g2, *fmt_args)
            drawing[g1y] += [chars.hbranch, label]
            positions[g2] = (g1x + 2, g1y)
            if g1x + 2 > width: width = g1x + 2
        else:
            # g1 has a descendant, so we branch g2 at the bottom of the drawing
            x, y = g1x + 2, len(drawing)
            label = fmt_node(g2, *fmt_args)
            drawing += [[chars.empty]*(g1x+1) + [chars.bottomleaf, label]]
            positions[g2] = (x,y)
            # modify characters above the branching to draw a nice branch.
            i = y - 1
            while i >= 0:
                if drawing[i][x-1] == chars.empty:
                    drawing[i][x-1] = chars.vbranch
                    i -= 1
                else:
                    if drawing[i][x-1] == chars.hbranch:
                        drawing[i][x-1] = chars.topnode
                    elif drawing[i][x-1] == chars.bottomleaf:
                        drawing[i][x-1] = chars.leaf
                    elif drawing[i][x-1] == chars.bottomnode:  # not possible
                        drawing[i][x-1] = chars.node
                    i = -1
    if '\n' in rootlabel:
        # split cells on multiple lines if needed
        for i in range(len(drawing)):
            line = drawing.pop(0)
            for newline in zip_longest(*(cell.split('\n') for cell in line),
                                       fillvalue=""):
                # need list to be able to assign in the next block
                drawing.append(list(newline))

        # optional but nicer: connect vertical branches
        for i, line in enumerate(drawing[:-1]):
            for j, cell in enumerate(line):
                if cell in [chars.topnode, chars.leaf, chars.vbranch] \
                        and drawing[i+1][j] == '':
                    drawing[i+1][j] = chars.vbranch

    # measure maximum width of each column
    colwidths = list(map(len, drawing[0]))
    for line in drawing[1:]:
        col = 0
        lcolwidths = len(colwidths)
        if lcolwidths < len(line):
            colwidths += map(len, line[lcolwidths:])
        while col < min(len(line), lcolwidths):
            w = len(line[col])
            if w > colwidths[col]:
                colwidths[col] = w
            col += 2
    
    out = ''
    regex_node = re.compile('[a-zA-Z0-9]')
    for line in drawing:
        # equalize column width
        for i, word in enumerate(line):
            fillchar = ' '
            if regex_node.search(word):
                fillchar = u'\u2500'
                #word = word.rstrip(' ')
            addchars = (colwidths[i]-len(word))
            if addchars:
                out += word + ' ' + fillchar * (addchars-1)
                if i+1 < len(line):
                    line[i+1] = fillchar + line[i+1].lstrip(' ')
            else:
                out += word
        out = out.rstrip(fillchar)
        out += '\n'
    return out.rstrip()

def draw_gene_trees(gene_trees, root, gene_info):
    fmt_node = lambda node, info: "%s (%.1s)" % (node, info[node].upper())
    return draw_trees(gene_trees, root, None, fmt_node, gene_info)


def draw_prot_trees(prottree):
    get_children = lambda treedata, node: [ch[0] for ch in treedata.get(node)]
    fmt_node = lambda node, info:\
                "%(taxon_name)-6.6s %(family_name)s (%(event).1s)" % info[node]
                #"%(family_name)s (%(event).1s)\n%(taxon_name)s" % info[node]
    return draw_trees(prottree.data, prottree.root, get_children, fmt_node,
                      prottree.info)


def draw_phyltree(phyltree):
    get_children = lambda phtr_items,node: [ch[0] for ch in \
                                            phtr_items.get(node, [])]
    fmt_node = lambda node: node
    return draw_trees(phyltree.items, phyltree.root, get_children, fmt_node)


def draw_nesteddict_tree(dicttree):
    #raise NotImplementedError
    get_children = lambda dicttree, node: list(node.keys())
    root = list(dicttree.keys())[0]
    return draw_trees(dicttree, root, get_children)


def draw_flatdict_tree(dicttree, root):
    return draw_trees(dicttree, root)


def draw_ete3_tree(ete3tree):
    get_children = lambda tree, node: node.children
    fmt_node = lambda node: node.name
    return draw_trees(ete3tree, ete3tree, get_children, fmt_node)  # *fmt_args)


def draw_newick(nwkfile, format=1):
    import ete3
    tree = ete3.Tree(nwkfile, format=format)
    tree.ladderize(direction=1)
    print(draw_ete3_tree(tree))


# Unfinished?
def draw_genetree_nice(gene_trees, root):
    positions = {}
    get_children_func = lambda tree,node: tree.get(node, set())
    lasty = [0] # current biggest y value per column (x value)
    drawing = [[],[]]
    for g, children in reversed(list(dfw_descendants_generalized(gene_trees,
                                        get_children_func, queue=[root]))):
        children_X = []
        children_Y = []
        for child in children:
            child_x, child_y = positions.get(child, (None, None))
            if not child_x and not child_y:
                # it is a terminal leaf
                child_x, child_y = -1, lasty[-1]
                lasty[-1] += 1
                drawing[child_x].append(chars.topleaf + "%-20s" % child)
            children_X.append(child_x)
            children_Y.append(child_y)
        children_y1, children_y2 = min(children_Y), max(children_Y)
        parent_Y = (children_y2 + children_y1) / 2
        parent_X = min(children_X) - 1
        if -len(drawing) > parent_X:
            drawing = [[], []] + drawing
        nodechar = chars.node
        if parent_Y == children_y1:
            if parent_Y == children_y2:
                nodechar = chars.hbranch
            else:
                nodechar = chars.topnode
        else:
            i = children_y1
            # fill vertically to the right start
            curr_y = len(drawing[parent_X])
            drawing[parent_X] += [chars.empty] * (i - curr_y)
            drawing[parent_X-1] += [chars.empty] * (i - curr_y)
            drawing[parent_X].append(chars.topleaf)
            drawing[parent_X - 1].append(chars.empty)
            while i < parent_Y:
                drawing[parent_X].append()

        for i in parent_Y: drawing[parent_X]


def make_lines(i, modgene, g, gene_categories, gene_trees, chrom=None, pos=1):
    cat = gene_categories.get(g)
    # add modern gene name as descendant of oldgene name
    ancmodgene = gene_trees.setdefault(g, set())
    modgene_i = modgene % i
    ancmodgene.add(modgene_i)

    # one chr per gene except tandems
    if chrom is None:
        chrom = modgene_i.replace('gene', 'chr')

    # lines for ANC, GO, GO_CAT, GEN
    ANC_line    = g + ' ' + modgene_i + '\n' # space as separator
    GO_line     = modgene_i + '\t' + str(cat) + '\n'
    GO_CAT_line = g + '\t' + GO_line
    GEN_line    = "%s\t%d\t%d\t1\t%s\t\n" % (chrom, pos, pos+1, modgene_i)
    return cat, chrom, (ANC_line, GO_line, GO_CAT_line, GEN_line)


def make_all_lines(genome, gene_categories, gene_trees, modgene, gene_info):
    """create ANC lines, GEN lines, OUT_GO lines, OUT_GO_CAT lines.
    Re-use the same chromosome name for adjacent genes
    Update genome, gene_trees"""
    i = 0

    while genome:
        g, neighbours = genome.popitem()
    
        i += 1
        cat, chrom, lines = make_lines(i, modgene, g, gene_categories,
                                       gene_trees)
        yield cat, lines

        # repeat again for adjacent genes
        pos = 1
        for ne in neighbours:
            try:
                genome.pop(ne)
            except KeyError as e:
                print(" -- ERROR -- ", file=sys.stderr)
                print("info:", gene_info[g], "g =", g, "n =", ne, file=sys.stderr)
                print("neighbors", neighbours, file=sys.stderr)
                gtree = g.split('.')[0]
                print(draw_gene_trees(gene_trees,gtree,gene_info), file=sys.stderr)
                raise e
            pos += 2
            i += 1
            cat, _, lines = make_lines(i, modgene, ne, gene_categories,
                                       gene_trees, chrom, pos)
            yield cat, lines

