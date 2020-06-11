#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from seqtools.plot_al_conservation import *
from Bio.Phylo.BaseTree import Clade
logging.basicConfig()
logger.setLevel(logging.DEBUG)

class Test_parsimony_between_parts:
    
    def test_score_dimension_is_1_and_shape_is_sequence_length(self):
        # sequences in rows. 0 is special: 
        alint = np.array([[1],[1],[2]])
        seqlabels = ['a', 'b', 'c']  # Where (a,b) forms a clade
        tree = Clade(name='r', clades=[Clade(name='ab',
                                             clades=[Clade(name='a'),
                                                     Clade(name='b')]),
                                       Clade(name='c')])
        score = Parsimony(alint, tree, seqlabels, minlength=4).rootwards()
        assert len(score.shape) == 1
        assert score.shape == (1,)

    def test_total_score_is_1_third(self):
        # sequences in rows. 0 is special: 
        alint = np.array([[1],[1],[2]])
        seqlabels = ['a', 'b', 'c']  # Where (a,b) forms a clade
        tree = Clade(name='r', clades=[Clade(name='ab',
                                             clades=[Clade(name='a'),
                                                     Clade(name='b')]),
                                       Clade(name='c')])
        score = Parsimony(alint, tree, seqlabels, minlength=4).rootwards()
        # There are 3 branches in the **unrooted** tree, and 1 substitution.
        assert np.isclose(score[0], 1./3)

    def test_part_score_returns_2_parts(self):
        # sequences in rows. 0 is special: 
        alint = np.array([[1],[1],[2]])
        seqlabels = ['a', 'b', 'c']  # Where (a,b) forms a clade
        tree = Clade(name='r', clades=[Clade(name='ab',
                                             clades=[Clade(name='a'),
                                                     Clade(name='b')]),
                                       Clade(name='c')])
        part_scores, part_branch_nbs = Parsimony(alint, tree, seqlabels, minlength=4,
                                                 parts=[(0,1)]).rootwards()
        assert len(part_scores) == 2
        assert len(part_branch_nbs) == 2

    def test_part_score_adds_stem_changes_to_outgroup(self):
        """Realistically, this behavior should be changed: changes should be
        added to the stem branch of the part clades.
        However this requires setting ancestral states, so an additional tree traversal."""
        # sequences in rows. 0 is special: 
        alint = np.array([[1],[1],[2]])
        seqlabels = ['a', 'b', 'c']  # Where (a,b) forms a clade
        tree = Clade(name='r', clades=[Clade(name='ab',
                                             clades=[Clade(name='a'),
                                                     Clade(name='b')]),
                                       Clade(name='c')])
        part_scores, part_branch_nbs = Parsimony(alint, tree, seqlabels, minlength=4,
                                              parts=[(0,1)]).rootwards()
        assert np.isclose(part_scores[0][0], 1)
        assert np.isclose(part_scores[1][0], 0)

    #def test_which_branch_gets_scored_at_the_root
        # The score should be "unrooted"

    def test_anc_states(self):
        alint = np.array([[1],[1],[2]])
        seqlabels = ['a', 'b', 'c']  # Where (a,b) forms a clade
        tree = Clade(name='r', clades=[Clade(name='ab',
                                             clades=[Clade(name='a'),
                                                     Clade(name='b')]),
                                       Clade(name='c')])
        parsimony = Parsimony(alint, tree, seqlabels, minlength=4)
        score_leafward = parsimony()
        print(parsimony.anc_states)
        assert not parsimony.anc_states[tree][0,0]  # not a gap
        assert parsimony.anc_states[tree][1,0]      # allowed to be state 1
        assert parsimony.anc_states[tree][2,0]      # allowed to be state 2
        clade_ab = tree.clades[0]
        assert not parsimony.anc_states[clade_ab][0,0]
        assert parsimony.anc_states[clade_ab][1,0]
        assert not parsimony.anc_states[clade_ab][2,0]  # not state 2
        
    def test_total_score_leafward(self):
        pass

