#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import itertools as it
from datasci.graphs import get_overlap_lims

# - a<b<c<d  no overlap
# - a<c<b<d  overlap
# - a<c<d<b  overlap (contained)
# +symetrics when range2 starts before:
# - c<d<a<b  no
# - c<a<d<b  overlap
# - c<a<b<d  overlap (contained)


class forOnePair(object):
    def test_number_of_ranges_is_1(self):
        (lims1, lims2, overlap) = get_overlap_lims(*self.get_ranges())
        print(lims1)
        print(lims2)
        print(overlap)
        assert lims1.shape == (2, 1)
        assert lims2.shape == (2, 1)
        assert overlap.shape == (2, 1)
    @classmethod
    def get_output(cls):
        #lims1, lims2, overlap = get_overlap_lims(*cls.shape_ranges(*cls.get_abcd()))
        lims1, lims2, overlap = get_overlap_lims(*cls.get_ranges())
        return lims1[:,0], lims2[:,0], overlap[:,0]

    # Should I do this...?
    @classmethod
    def get_abcd(cls):
        letters = list(cls.__name__.rsplit('_', 1)[-1])
        return [1+letters.index(l) for l in 'abcd']

    @staticmethod
    def shape_ranges(a,b,c,d):
        return [[a],[b]], [[c],[d]]

class for3Pairs(object):
    def test_input_ranges_are_3x2(self):
        # Testing the test.
        ranges_ab, ranges_cd = self.get_ranges()
        assert np.asarray(ranges_ab).shape == (2,3)
        assert np.asarray(ranges_cd).shape == (2,3)

    def test_number_of_ranges_is_3(self):
        (lims1, lims2, overlap) = get_overlap_lims(*self.get_ranges())
        print(lims1)
        print(lims2)
        print(overlap)
        assert lims1.shape == (2, 3)
        assert lims2.shape == (2, 3)
        assert overlap.shape == (2, 3)
    @classmethod
    def get_output(cls):
        #lims1, lims2, overlap = get_overlap_lims(*cls.shape_ranges(*cls.get_abcd()))
        return get_overlap_lims(*cls.get_ranges())

    # Should I do this...?
    @classmethod
    def get_abcd(cls):
        letters = list(cls.__name__.rsplit('_', 1)[-1])
        return [1+letters.index(l) for l in 'abcd']

    @staticmethod
    def shape_ranges(a,b,c,d):
        return [[a]*3,[b]*3], [[c]*3,[d]*3]
    

class noOverlap(object):
    def test_overlap_is_nan(self):
        (lims1, lims2, overlap) = self.get_output()
        assert np.isnan(overlap).all()

    def test_lim1_unchanged(self):
        range_ab, range_cd = self.get_ranges()
        (lims1, lims2, overlap) = self.get_output()
        assert lims1.astype(int).tolist() == np.array(range_ab).flatten().astype(int).tolist()

    def test_lim2_unchanged(self):
        range_ab, range_cd = self.get_ranges()
        (lims1, lims2, overlap) = self.get_output()
        print(range_ab, range_cd, lims1, lims2, overlap, sep='\n')
        assert lims2.astype(int).tolist() == np.array(range_cd).flatten().astype(int).tolist()

class startsWith1(object):
    def test_lim1_startswith_1(self):
        (lims1, lims2, overlap) = self.get_output()
        assert lims1[0] == 1

class startsWith2(object):
    def test_lim2_startswith_1(self):
        (lims1, lims2, overlap) = self.get_output()
        assert lims2[0] == 1


class Test_abcd(forOnePair, noOverlap, startsWith1):
    @classmethod
    def get_ranges(cls):
        a, b, c, d = 1, 2, 3, 4
        return cls.shape_ranges(a,b,c,d)


class Test_cdab(forOnePair, noOverlap, startsWith2):
    @classmethod
    def get_ranges(cls):
        c, d, a, b = 1, 2, 3, 4
        return cls.shape_ranges(a,b,c,d)


class partialOverlap(object):
    pass

class containedOverlap(object):
    #@abstractmethod
    #def test_contained_is_nan(self):
    #def test_overlap_is_cd_or_ab
    pass

class Test_acbd(forOnePair,partialOverlap,startsWith1):
    @classmethod
    def get_ranges(cls):
        a, c, b, d = 1, 2, 3, 4
        return cls.shape_ranges(a,b,c,d)
    def test_lim1_is_ac(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert lims1.astype(int).tolist() == [a, c]
    def test_overlap_is_cb(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert overlap.astype(int).tolist() == [c, b]
    def test_lim2_is_bd(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert lims2.astype(int).tolist() == [b, d]

class Test_cadb(forOnePair,partialOverlap,startsWith2):
    @classmethod
    def get_ranges(cls):
        c, a, d, b = 1, 2, 3, 4
        return cls.shape_ranges(a,b,c,d)
    def test_lim1_is_db(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert lims1.astype(int).tolist() == [d, b]
    def test_overlap_is_ad(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert overlap.astype(int).tolist() == [a, d]
    def test_lim2_is_ca(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert lims2.astype(int).tolist() == [c, a]

class Test_acdb(forOnePair,containedOverlap,startsWith1):
    @classmethod
    def get_ranges(cls):
        a, c, d, b = 1, 2, 3, 4
        return cls.shape_ranges(a,b,c,d)
    def test_lim1_is_ab(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert lims1.astype(int).tolist() == [a, b]
    def test_lim2_is_nan(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert np.isnan(lims2).all()
    def test_overlap_is_cd(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert overlap.astype(int).tolist() == [c, d]


class Test_cabd(forOnePair,containedOverlap,startsWith2):
    @classmethod
    def get_ranges(cls):
        c, a, b, d = 1, 2, 3, 4
        return cls.shape_ranges(a,b,c,d)
    def test_lim1_is_nan(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert np.isnan(lims1).all()
    def test_lim2_is_cd(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert lims2.astype(int).tolist() == [c, d]
    def test_overlap_is_ab(self):
        a, b, c, d = self.get_abcd()
        (lims1, lims2, overlap) = self.get_output()
        assert overlap.astype(int).tolist() == [a, b]


class Test_3pairs_abcd(for3Pairs):
    @classmethod
    def get_ranges(cls):
        return cls.shape_ranges(*cls.get_abcd())


class Test_6pairs_all_possible_combinations(object):
    @classmethod
    def get_ranges(cls):
        # This generate impossible combinations, with b<a or d<c
        ranges = np.array([p for p in it.permutations(range(4), 4)
                            if (p[0]<p[1] and p[2]<p[3])]).T
        return ranges[:2].copy(), ranges[2:].copy()
        
    def test_all_combinations_in_one_list(self):
        lims1, lims2, overlap = get_overlap_lims(*self.get_ranges())
