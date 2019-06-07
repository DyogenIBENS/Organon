#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from io import StringIO

from dendro.dollocorr import *


# NB: All asserted values have been manually calculated with pen & paper


# Test the child sum combination counters:
c1,c2,c3 = [1,1], [1,1], [1,1]
# 2 children first, and equality of the general and specific method.
assert 1 == combinations_of_2_children_summing_to_k([c1, c2], 0)
assert 1 == combinations_of_children_summing_to_k([c1, c2], 0)
assert 2 == combinations_of_2_children_summing_to_k([c1, c2], 1)
assert 2 == combinations_of_children_summing_to_k([c1, c2], 1)
assert 1 == combinations_of_2_children_summing_to_k([c1, c2], 2)
assert 1 == combinations_of_children_summing_to_k([c1, c2], 2)

assert 1 == combinations_of_children_summing_to_k([c1, c2, c3], 0)
assert 3 == combinations_of_children_summing_to_k([c1, c2, c3], 1)
assert 3 == combinations_of_children_summing_to_k([c1, c2, c3], 2)
assert 1 == combinations_of_children_summing_to_k([c1, c2, c3], 3)

c1 = [1,3]
assert 4 == combinations_of_2_children_summing_to_k([c1, c2], 1)
assert 4 == combinations_of_children_summing_to_k([c1, c2], 1)
assert 3 == combinations_of_2_children_summing_to_k([c1, c2], 2)
assert 3 == combinations_of_children_summing_to_k([c1, c2], 2)

assert 5 == combinations_of_children_summing_to_k([c1, c2, c3], 1)
assert 0+4+3 == combinations_of_children_summing_to_k([c1, c2, c3], 2)

c1, c2 = [1,3,1], [1,3,1]  # Possibilities for a tree of 2 leaves.
assert 6 == combinations_of_2_children_summing_to_k([c1, c2], 1)
assert 6 == combinations_of_children_summing_to_k([c1, c2], 1)
assert 11 == combinations_of_2_children_summing_to_k([c1, c2], 2)
assert 11 == combinations_of_children_summing_to_k([c1, c2], 2)

# Forbid some values:
c1, c2 = [0, 3, 1], [0, 3, 1]
assert 0 == combinations_of_2_children_summing_to_k([c1, c2], 0)
assert 0 == combinations_of_children_summing_to_k([c1, c2], 0)
assert 0 == combinations_of_2_children_summing_to_k([c1, c2], 1)
assert 0 == combinations_of_children_summing_to_k([c1, c2], 1)
assert 9 == combinations_of_2_children_summing_to_k([c1, c2], 2)
assert 9 == combinations_of_children_summing_to_k([c1, c2], 2)

# Smallest tree: one node.
tree1 = myPhylTree.PhylogeneticTree(StringIO("r;"))
r = place_single_events(1, tree1)

assert (r == [1,0]).all()

tree2 = myPhylTree.PhylogeneticTree(StringIO("(A:1,B:1)r;"))
r0 = place_single_events(0, tree2)

#assert isinstance(r0, dict)
assert r0[0] == 1
assert len(r0) <= 1

r1 = place_single_events(1, tree2)

assert (r1 == [1,2]).all()

r2 = place_single_events(2, tree2)

assert (r2 == [1,2,1]).all()

r3 = place_single_events(3, tree2)

assert (r3 == [1,2,1,0]).all()

# Smallest (interesting) tree: 3 leaves

tree3 = myPhylTree.PhylogeneticTree(StringIO("((A:1,B:1)ab:1,C:2)r;"))
r = place_single_events(1, tree3)

assert (r == [1,4]).all()

r = place_single_events(4, tree3)

assert (r == [1,4,4,1,0]).all()

tree4balanced = myPhylTree.PhylogeneticTree(StringIO("((A:1,B:1)ab:1,(C:0.5,D:0.5)cd:1.5)r;"))
r = place_single_events(6, tree4balanced)

assert (r == [1,6,11,6,1,0,0]).all()

tree4caterpillar = myPhylTree.PhylogeneticTree(StringIO("(((A:1,B:1)ab:0.5,C:1.5)abc:0.5,D:2.5)r;"))
r = place_single_events(6, tree4caterpillar)

assert (r == [1,6,9,5,1,0,0]).all()

# A node with 1 child.
# r -- a -- A
#   \- B
tree_nondicho = myPhylTree.PhylogeneticTree(StringIO("((A:1)a:1,B:2)r;"))
r = place_single_events(3, tree_nondicho)

assert (r == [1, 3, 2, 0]).all()

rb = place_single_events(3, myPhylTree.PhylogeneticTree(StringIO("(A:1,B:2)r;")))

# A node with 3 children.
tree_tricho = myPhylTree.PhylogeneticTree(StringIO("(A:1,B:2,C:1.5)r;"))
r = place_single_events(3, tree_tricho)

assert (r == [1,3,3,1]).all()


p1more, ptot = maddison_test(1, 2, tree4balanced, ['ab'])

assert p1more == 10
assert ptot == 11

p1more, ptot = maddison_test(1, 2, tree4caterpillar, ['ab'])

assert p1more == 7
assert ptot == 9

# Fig. 2 from Maddison 1990 (Correlated evolution of two binary characters)
tree_Maddison2 = myPhylTree.PhylogeneticTree(StringIO("((((a,b)ab,c)abc,d)G,((u,v)uv,(w,x)wx)H)F;"))
# From the article: 69 ways to have 2 gains (and zero loss)
assert 69 == place_single_events(2, tree_Maddison2)[2]


