#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import argparse
from scipy.stats import gamma


cladename_re = re.compile(r'<Gamma id="Gamma\.([A-Za-z_]+)" .*offset="([0-9.]+)"')
alpha_re = re.compile(r'<parameter id="RealParameter\.[0-9]+" .*name="alpha">([0-9.]+)</parameter>')
beta_re = re.compile(r'<parameter id="RealParameter\.[0-9]+" .*name="beta">([0-9.]+)</parameter>')


def check_gamma_param(paramfile, taxonagesfile):
    params = {}
    
    previous = [None]*3
    with open(paramfile) as f:
        for line in f:
            m = cladename_re.search(line)
            if m:
                clade, offset = m.groups()
                loc = float(offset)
                previous[0] = loc
            else:
                m_alpha = alpha_re.search(line)
                if m_alpha:
                    a = float(m_alpha.group(1))
                    previous[1] = a
                m_beta = beta_re.search(line)
                if m_beta:
                    scale = float(m_beta.group(1))
                    previous[2] = scale
            if all(p is not None for p in previous):
                params[clade] = previous
                previous = [None]*3

    with open(taxonagesfile) as f:
        for line in f:
            taxon, age = line.rstrip().split('\t')
            age = float(age)

            for clade, (loc, a, scale) in params.items():
                if taxon.startswith(clade):
                    print('%-35s -> %-20s    ' %(taxon, clade) +
                          '%5.2f + G(%5.2f,%5.2f) P(x < %5.2f) =' % (loc, a, scale, age),
                            gamma.cdf(age, a, loc, scale))
                    break
            else:
                print('Taxon %r NOT found.' % taxon)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paramfile')
    parser.add_argument('taxonagesfile')
    
    args = parser.parse_args()
    check_gamma_param(**vars(args))

