#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from codeml.prune2family import *
from dendro.reconciled import infer_gene_event


plogger = logger
plogger.setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


ensembl_version = 93


sample_data = [
        {'description': 'Regular speciation',
         'rootname': 'PanENSGT00760000119060.a',
         'treedata': '(ENSPTRG0:0.05,ENSPPAG0:0.04){0}:0.1;',
         'checks': {'insert_species_nodes_back': dict(nleaves=2,
                       ndescendants=2,
                       nreinserted=0,
                       childset=set(('ENSPTRG0', 'ENSPPAG0')),
                       event='spe')}
        },
        {'description': 'Missing species -> missing speciation',
         'rootname': 'HomoPanENSGT00760000119060.a',
         'treedata': '(ENSG000:0.05,ENSPPAG0:0.04){0}:0.1;',
         'checks': {'insert_species_nodes_back': dict(
                        nleaves=2,
                        ndescendants=3,
                        nreinserted=1,
                        childset=set(('ENSG000', 'PanENSGT00760000119060.a')),
                        event='spe')}
        },
        {'description': 'Regular duplication.',
         'rootname': 'Pan.paniscusENSGT00760000119060.a',
         'treedata': '(ENSPPAG0:0.05,ENSPPAG1:0.04){0}:0.1;',
         'checks': {'insert_species_nodes_back': dict(
                        nleaves=2,
                        ndescendants=2,
                        nreinserted=0,
                        childset=set(('ENSPPAG0', 'ENSPPAG1')),
                        event='dup')}
        },
        {'description': 'Missing species -> missing dup + speciation',
         'rootname': 'PanENSGT00760000119060.a',
         'treedata': '((ENSPPAG0:0.01,ENSPTRG0:0.01){0}.a:0.05,ENSPPAG1:0.04){0}:0.1;',
         'checks': {'insert_species_nodes_back': dict(
                        nleaves=3,
                        ndescendants=5,
                        nreinserted=1,
                        childset=set(('PanENSGT00760000119060.a.a',
                                      'PanENSGT00760000119060.a.b')),
                        event='implicit dup+spe')}
        },
        {'description': 'Not direct duplication (many missing speciation)',
         'rootname': 'HomoPanENSGT00760000119060.a',
         'treedata': '(ENSPPAG0:0.01,ENSPPAG1:0.04){0}:0.1;',
         'checks': {'insert_species_nodes_back': dict(
                        nleaves=2,
                        ndescendants=6,
                        nreinserted=4,
                        #childset=set(()),
                        event='implicit dup+spe')}
        },
        {'description': 'Not direct duplication (one missing speciation)',
         'rootname': 'PanENSGT00760000119060.a',
         'treedata': '(ENSPPAG0:0.01,ENSPPAG1:0.04){0}:0.1;',
         'checks': {'insert_species_nodes_back': dict(
                        nleaves=2,
                        ndescendants=4,
                        nreinserted=2,
                        #childset=set(()),
                        event='implicit dup+spe')}
        },
        {'description': 'Not direct speciation (many missing speciation)',
         'rootname': 'HomoPanENSGT00760000119060.a',
         'treedata': '(ENSPTRG0:0.01,ENSPPAG0:0.04){0}:0.1;',
         'checks': {'insert_species_nodes_back': dict(
                        nleaves=2,
                        ndescendants=4,
                        nreinserted=2,
                        #childset=set(()),
                        event='spe not mrca')}
        },
        # 7.
        {'description': 'Distant outgroup',
         'rootname': 'HomoPanENSGT00760000119060.a',
         'treedata': '((ENSPPAG0:0.01,ENSPPAG1:0.04)Pan.paniscusENSGT00760000119060:0.1,ENSG000:0.3){0}:0.1;',
         'checks': {'insert_species_nodes_back': dict(
                        nleaves=3,
                        ndescendants=5,
                        nreinserted=1,
                        #childset=set(()),
                        event='spe')}
        }]


# Non treebest mode:
ancgene2sp = re.compile(r'(Pan\.troglodytes|Pan\.paniscus|Homo\.sapiens|HomoPan|Pan)([^a-z].*|)$')
def get_species(node):
    return ultimate_seq2sp(node.name, ensembl_version), node.name
def split_ancestor(node):
    return split_species_gene(node.name, ancgene2sp)


this_parse_species_genename = partial(parse_species_genename,
                                      get_species=get_species,
                                      split_ancestor=split_ancestor)

diclinks = {'HomoPan': {'Homo sapiens': ['HomoPan', 'Homo sapiens'],
                        'Pan': ['HomoPan', 'Pan'],
                        'Pan paniscus': ['HomoPan', 'Pan', 'Pan paniscus'],
                        'Pan troglodytes': ['HomoPan', 'Pan', 'Pan troglodytes'],
                        'HomoPan': ['HomoPan']},
            'Pan': {'Pan paniscus': ['Pan', 'Pan paniscus'],
                    'Pan troglodytes': ['Pan', 'Pan troglodytes'],
                    'Pan': ['Pan']},
            'Homo sapiens': {'Homo sapiens': ['Homo sapiens']},
            'Pan troglodytes': {'Pan troglodytes': ['Pan troglodytes']},
            'Pan paniscus': {'Pan paniscus': ['Pan paniscus']}
            }

ages = {'HomoPan': 6.65, 'Pan': 2.82, 'Homo sapiens': 0, 'Pan paniscus': 0,
        'Pan troglodytes': 0}


def test_insert_species_nodes_back(
                                   fix_suffix=True,
                                   force_mrca=False):
    def run(treedata):
        tree = ete3.Tree(treedata, format=1)
        children_taxa = []
        for child in tree.children:
            try:
                children_taxa.append(ultimate_seq2sp(child.name, ensembl_version))
            except KeyError:
                children_taxa.append(split_ancestor(child)[0])

        print('infer_gene_event:',
              infer_gene_event(tree,
                              split_ancestor(tree)[0],
                              set(children_taxa)))

        insert_species_nodes_back(tree, this_parse_species_genename, diclinks, ages,
                                  fix_suffix, force_mrca)
        return tree

    def check(tree, rootname, nleaves, ndescendants, nreinserted, event,
              childset=None):
        try:
            assert tree.name == rootname, tree.name
            assert len(tree) == nleaves, len(tree)
            assert len(tree.get_descendants()) == ndescendants, len(tree.get_descendants())
            assert sum(int(getattr(n, 'reinserted', None) is not None)
                        for n in tree.traverse()) == nreinserted
            assert tree.event == event, tree.event
            if childset is not None:
                assert set(ch.name for ch in tree.children) == childset
        except AssertionError as err:
            logger.exception('Failed check:')
            print('occured on:\n' + tree.get_ascii())
            return False
        return True

    passed = []
    total = 0

    for sample in sample_data:
        check_values = sample['checks'].get('insert_species_nodes_back')
        if check_values is not None:
            total += 1
            print('# %d. %s.' % (total, sample['description']))
            tree = run(sample['treedata'].format(sample['rootname']))
            result = check(tree, sample['rootname'], **check_values)
            print(result)
            passed.append(int(result))

    print('-' * 20)
    print('%s:passed %d/%d tests: [%s]' % (
            'insert_species_nodes_back', sum(passed), total,
            ' '.join(str(p) for p in passed)))
    print('-' * 20)
    return sum(passed)==total


def test_reroot_with_outgroup():

    global print_if_verbose
    def print_if_verbose(*args, **kwargs):
        pass

    def check(orig_tree, tree, outsize, treetype, rootname, expect_outsize, outgroupnodes, childset=None):
        try:
            assert not orig_tree.search_nodes(is_outgroup=1), "Left outgroup marks."
            assert not orig_tree.search_nodes(has_ingroup=1), "Left ingroup marks."
            if treetype is None:
                assert tree is None, type(tree)
            else:
                assert isinstance(tree, ete3.TreeNode)

                assert tree.name == rootname, tree.name
                marked_outgroups = set(n.name for n in tree.search_nodes(is_outgroup=1))
                assert marked_outgroups == outgroupnodes, 'marked %s VS predicted %s' %(marked_outgroups, outgroupnodes)
                if childset is not None:
                    assert set(ch.name for ch in root.children) == childset

            assert outsize == expect_outsize, outsize
        except AssertionError as err:
            logger.exception('Failed check:')
            print('occured on:\n', 'None' if tree is None else tree.write(features=['is_outgroup', 'has_ingroup', 'event'], format=1, format_root_node=True))
            return False
        return True

    sample = sample_data[7]
    tree = ete3.Tree(sample['treedata'].format(sample['rootname']), format=1)
    insert_species_nodes_back(tree, this_parse_species_genename, diclinks, ages,
                              fix_suffix=True, force_mrca=False)
    
    passed = []
    total = 0

    print(tree.get_ascii())
    for leaf in tree.iter_leaves():
        total += 1
        print('> TREE:', tree.write(features=['is_outgroup', 'has_ingroup', 'event'], format=1, format_root_node=True))
        
        # Predicted result  (Works because the tree has 3 leaves)
        outgroupnodes = [l.name for l in tree.iter_leaves() if l != leaf]
        common_outgroup = tree.get_common_ancestor(outgroupnodes)
        if common_outgroup != tree:
            outgroupnodes = [common_outgroup.name]

        root, outsize = reroot_with_outgroup(leaf, maxsize=2, minsize=2,
                                             is_allowed_outgroup=None)
        if root is not None:
            print('> old node %r (%s) VS new node %r (%s)' % (leaf.name, id(leaf),
                            root.name, id(root)))
        result = check(tree, root, outsize, treetype=ete3.TreeNode,
                       rootname='HomoPanENSGT00760000119060.a', # tree.name
                       expect_outsize=2,
                       outgroupnodes=set(outgroupnodes))

        print('Leaf test 1:', result)
        passed.append(int(result))

        total += 1
        root, outsize = reroot_with_outgroup(leaf, maxsize=3, minsize=3,
                                             is_allowed_outgroup=None)
        if root is not None:
            print('> old node %r (%s) VS new node %r (%s)' % (leaf.name, id(leaf),
                            root.name, id(root)))
        result = check(tree, root, outsize, treetype=None,
                       rootname=None, # tree.name
                       expect_outsize=2,
                       outgroupnodes=set(outgroupnodes))

        print('Leaf test 2:', result)
        passed.append(int(result))

    ##
    tree = ete3.Tree('(a,b)r;', format=1)
    root, outsize = reroot_with_outgroup(tree&'a', maxsize=2, minsize=2,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=None, rootname=None,
                    expect_outsize=1, outgroupnodes=set())
    total += 1
    passed.append(int(result))

    ##
    root, outsize = reroot_with_outgroup(tree&'a', maxsize=1, minsize=1,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=ete3.TreeNode, rootname='r',
                    expect_outsize=1, outgroupnodes=set('b'))
    total += 1
    passed.append(int(result))

    ## Check that the original tree was not altered by the call, and that the outgroup of another node can be found.
    root, outsize = reroot_with_outgroup(tree&'b', maxsize=1, minsize=1,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=ete3.TreeNode, rootname='r',
                    expect_outsize=1, outgroupnodes=set('a'))
    total += 1
    passed.append(int(result))

    ## Add a single-child node.
    #  minsize=1 and minsize=0 should behave differently -> answer NO.
    #  But minsize=2 and minsize=1 should.
    tree = ete3.Tree('((a)aa,b)r;', format=1)
    root, outsize = reroot_with_outgroup(tree&'a', maxsize=1, minsize=1,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=ete3.TreeNode, rootname='r',
                    expect_outsize=1, outgroupnodes=set('b'))
    total += 1
    passed.append(int(result))

    root, outsize = reroot_with_outgroup(tree&'a', maxsize=1, minsize=0,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=ete3.TreeNode, rootname='r',
                    expect_outsize=1, outgroupnodes=set('b'))
    total += 1
    passed.append(int(result))
    
    ##
    tree = ete3.Tree('(((a)aa,b)x,c)r;', format=1)
    root, outsize = reroot_with_outgroup(tree&'a', maxsize=1, minsize=1,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=ete3.TreeNode, rootname='x',
                    expect_outsize=1, outgroupnodes=set('b'))
    total += 1
    passed.append(int(result))
    ## Check from the node b:
    root, outsize = reroot_with_outgroup(tree&'b', maxsize=1, minsize=1,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=ete3.TreeNode, rootname='x',
                    expect_outsize=1, outgroupnodes=set('a'))
    total += 1
    passed.append(int(result))
    root, outsize = reroot_with_outgroup(tree&'b', maxsize=2, minsize=2,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=ete3.TreeNode, rootname='r',
                   expect_outsize=2, outgroupnodes=set('ac'))
    total += 1
    passed.append(int(result))

    # Test the thinning.
    root, outsize = reroot_with_outgroup(tree&'c', maxsize=1, minsize=1,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=ete3.TreeNode, rootname='r',
                   expect_outsize=1, outgroupnodes=set('b'))
    total += 1
    passed.append(int(result))

    tree = ete3.Tree('(((a)aa,b)x,c)r;', format=1)
    root, outsize = reroot_with_outgroup(tree&'a', maxsize=2, minsize=1,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=ete3.TreeNode, rootname='x',
                    expect_outsize=1, outgroupnodes=set('b'))
    total += 1
    passed.append(int(result))

    root, outsize = reroot_with_outgroup(tree&'a', maxsize=2, minsize=2,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=ete3.TreeNode, rootname='r',
                    expect_outsize=2, outgroupnodes=set('bc'))
    total += 1
    passed.append(int(result))

    ## One more single-child node.
    tree = ete3.Tree('((((a)aa,b)x)xx,c)r;', format=1)
    root, outsize = reroot_with_outgroup(tree&'a', maxsize=2, minsize=2,
                                         is_allowed_outgroup=None)
    result = check(tree, root, outsize, treetype=ete3.TreeNode, rootname='r',
                    expect_outsize=2, outgroupnodes=set('bc'))
    total += 1
    passed.append(int(result))

    ## END: summary
    print('-' * 20)
    print('%s:passed %d/%d tests: [%s]' % ('reroot_with_outgroup',
            sum(passed), total, ' '.join(str(p) for p in passed)))
    print('-' * 20)
    return sum(passed)==total


if __name__=='__main__':
    logging.basicConfig(format=logging.BASIC_FORMAT)

    r = True
    r &= test_insert_species_nodes_back()
    r &= test_reroot_with_outgroup()
    
    exit(0 if r else 1)

