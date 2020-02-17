#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import ete3
from dendro.any import skip_set_children, ete3 as ete3_methods
import dendro.sorter as ds

get_children = ete3_methods.get_children
set_children = ete3_methods.set_children
def verbose_set_children(tree, node, children):
    print('Assigning to %s -> %s' % (node.name, [ch.name for ch in children]))
    set_children(tree, node, children)

class Test_assign_children(object):
    def test_tree_a(self):
        tree = ete3.Tree('(a,b)r;', format=1)
        old_children = tree.children
        assert [ch.name for ch in old_children] == ['a', 'b']
        verbose_set_children(tree, tree, list(reversed(old_children)))
        new_children = tree.children
        assert [ch.name for ch in new_children] == ['b', 'a']
    def test_tree_a_getchildren_method(self):
        tree = ete3.Tree('(a,b)r;', format=1)
        old_children = tree.get_children()
        assert [ch.name for ch in old_children] == ['a', 'b']
        verbose_set_children(tree, tree, list(reversed(old_children)))
        new_children = tree.get_children()
        assert [ch.name for ch in new_children] == ['b', 'a']
    def test_tree_a_getchildren_func(self):
        tree = ete3.Tree('(a,b)r;', format=1)
        old_children = get_children(tree, tree)
        assert [ch.name for ch in old_children] == ['a', 'b']
        verbose_set_children(tree, tree, list(reversed(old_children)))
        new_children = get_children(tree, tree)
        assert [ch.name for ch in new_children] == ['b', 'a']
    def test_tree_a_assign_same_children_copy(self):
        tree = ete3.Tree('(a,b)r;', format=1)
        old_children = tree.children
        assert [ch.name for ch in old_children] == ['a', 'b']
        verbose_set_children(tree, tree, list(old_children))
        new_children = tree.children
        assert [ch.name for ch in new_children] == ['a', 'b']
        assert old_children == new_children
    def test_tree_a_assign_same_children(self):
        tree = ete3.Tree('(a,b)r;', format=1)
        old_children = tree.children
        assert [ch.name for ch in old_children] == ['a', 'b']
        verbose_set_children(tree, tree, old_children)
        new_children = tree.children
        assert [ch.name for ch in new_children] == ['a', 'b']
        assert old_children == new_children


class Test_cumul_mass(object):
    """Check that it returns the number of leaves below a node."""
    def test_tree_a(self):
        tree = ete3.Tree('(a,b)r;', format=1)
        sizes = ds.ladderize(tree, tree, get_children, assign=skip_set_children)
        assert sizes[tree&'a'] == 1
        assert sizes[tree&'b'] == 1
        assert sizes[tree&'r'] == 2
    def test_tree_c(self):
        tree = ete3.Tree('(((c1, c2)c,d)a,b)r;', format=1)
        sizes = ds.ladderize(tree, tree, get_children, assign=skip_set_children)
        assert sizes[tree&'b'] == 1
        assert sizes[tree&'c1'] == 1
        assert sizes[tree&'c2'] == 1
        assert sizes[tree&'c'] == 2
        assert sizes[tree&'a'] == 3
        assert sizes[tree&'r'] == 4


class Test_find_pivot(object):
    
    def test_tree_a(self):
        tree = ete3.Tree('(a,b)r;', format=1)
        pivot, split_index = ds.find_pivot(tree, tree, get_children, set_children)
        assert pivot == tree
        assert split_index == 1
    def test_tree_b(self):
        tree = ete3.Tree('((a1,a2)a,b)r;', format=1)
        pivot, split_index = ds.find_pivot(tree, tree, get_children, verbose_set_children)
        assert pivot == tree, 'Should have stopped at node "r"'
        assert split_index == 1
        assert tree.children[0].name == 'b', "'b' should be the top child of 'r'"
    def test_tree_c(self):
        tree = ete3.Tree('(((c1, c2)c,d)a,b)r;', format=1)
        #         /-c1
        #      /c|
        #   /a|   \-c2
        #  |  |
        #-r|   \-d
        #  |
        #   \-b
        # Output should be:
        #    
        #    
        #   /-b
        #  | 
        #-r|   /-d
        #  |  |<<
        #   \a|   /-c1
        #      \c|
        #         \-c2
        #     
        
        pivot, split_index = ds.find_pivot(tree, tree, get_children, set_children)
        assert pivot == (tree&'a')
        assert split_index == 1
        assert (tree&'a').children[0].name == 'd'
    def test_tree_equal_deltas_choose_previous(self):
        tree = ete3.Tree('(((d1,d2)d,c)a,(b1,b2)b)r;', format=1)
        #          /-d1
        #       /d|
        #    /a|   \-d2
        #   |  |        <--
        # -r|   \-c
        #   |           <-- choose here.
        #   |   /-b1
        #    \b|
        #       \-b2

        # We could split at 'r', or at 'a' with the same result. The algo should choose 'r'.
        pivot, split_index = ds.find_pivot(tree, tree, get_children, set_children)
        assert pivot == tree
        assert split_index == 1
    def test_tree_previous_delta_the_smallest(self):
        #          /-1
        #       /d|--2
        #    /a|   \-3
        #   |  |        <--- 
        #   |  |   /-1
        #   |   \c|
        # -r|      \-2
        #   |           <--- choose here because 4 leaves ingroup VS 5 leaves outgroup.
        #   |   /-1
        #   |  |--2
        #    \b|--3
        #       \-4

        tree = ete3.Tree('(((1,2,3)d,(1,2)c)a,(1,2,3,4)b)r;', format=1)
        pivot, split_index = ds.find_pivot(tree, tree, get_children, set_children)
        assert pivot == tree
    def test_tree_next_delta_the_smallest(self):
        #         /-1
        #        |--2
        #      /d|--3
        #     |   \-4
        #   /a|        <--- choose here
        #  |  |   /-1
        #  |   \c|
        #-r|      \-2
        #  |           <---
        #  |   /-1
        #   \b|--2
        #      \-3

        tree = ete3.Tree('(((1,2,3,4)d,(1,2)c)a,(1,2,3)b)r;', format=1)
        pivot, split_index = ds.find_pivot(tree, tree, get_children, set_children)
        assert pivot == (tree&'a')


if __name__ == '__main__':
    from sys import exit
    results = []
    #print(globals())
    #print(vars(Test_find_pivot))
    for testclass in [obj for name,obj in globals().items() if name.startswith('Test')]:
        testcase = testclass()
        #print(dir(testcase))
        for testname in [attrname for attrname in dir(testcase) if
                         attrname.startswith('test_')]:
            try:
                getattr(testcase, testname)()
                results.append(1)
            except AssertionError:
                results.append(0)

    #testcase = Test_find_pivot()
    #testcase.test_tree_a()
    #testcase.test_tree_b()
    #testcase.test_tree_c()
    
    print(results)
    if not all(results):
        exit(1)

