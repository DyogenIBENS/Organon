#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from LibsDyogen import myPhylTree
from io import StringIO

from dollocorr import *

# Smallest tree: one node.
tree1 = myPhylTree.PhylogeneticTree(StringIO("r;"))
r = place_single_events(1, tree1)

assert isinstance(r, dict)
assert r[0] == 1
assert r[1] == 0 

tree2 = myPhylTree.PhylogeneticTree(StringIO("(A:1,B:1)r;"))
r0 = place_single_events(0, tree2)

#assert isinstance(r0, dict)
assert r0[0] == 1
assert 1 not in r0.keys()

r1 = place_single_events(1, tree2)

assert r1[0] == 1
assert r1[1] == 2

r2 = place_single_events(2, tree2)

assert r2[0] == 1
assert r2[1] == 2
assert r2[2] == 1

r3 = place_single_events(3, tree2)

assert r3[0] == 1
assert r3[1] == 2
assert r3[2] == 1
#assert r3[3] == 0

# Smallest (interesting) tree: 3 leaves

tree3 = myPhylTree.PhylogeneticTree(StringIO("((A:1,B:1)ab:1,C:2)r;"))
r = place_single_events(1, tree3)

assert r[0] == 1
assert r[1] == 4

r = place_single_events(4, tree3)

assert r[0] == 1
assert r[1] == 4
assert r[2] == 4
assert r[3] == 1
assert r[4] == 0

tree4balanced = myPhylTree.PhylogeneticTree(StringIO("((A:1,B:1)ab:1,(C:0.5,D:0.5)cd:1.5)r;"))
r = place_single_events(6, tree4balanced)

assert r[0] == 1
assert r[1] == 6
assert r[2] == 11
assert r[3] == 6
assert r[4] == 1
assert r[5] == 0
assert r[6] == 0

tree4caterpillar = myPhylTree.PhylogeneticTree(StringIO("(((A:1,B:1)ab:0.5,C:1.5)abc:0.5,D:2.5)r;"))
r = place_single_events(6, tree4caterpillar)

assert r[0] == 1
assert r[1] == 6
assert r[2] == 9
assert r[3] == 5
assert r[4] == 1
assert r[5] == 0
assert r[6] == 0

tree_nondicho = myPhylTree.PhylogeneticTree(StringIO("((A:1)a:1,B:2)r;"))
try:
    r = place_single_events(1, tree_nondicho)
except ValueError as expected_err:
    assert 'Not implemented for non dichotomic trees' in expected_err.args[0]
