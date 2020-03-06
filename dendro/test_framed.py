#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from dendro.framed import *


tree = ete3.Tree('(((e:0.5,f)c,d)a,(g:2,h:2)b)r;', format=1)
tree2 = ete3.Tree('(((e2:0.5,f2)d2,c2)a2,(g2:2,h2:2)b2)r2;', format=1)

treedf = ete3_to_parentdata(tree)
treedf2 = ete3_to_parentdata(tree2)

trees_df = pd.concat((treedf, treedf2))

leaves = ['e', 'f', 'd', 'g', 'h', 'e2', 'f2', 'c2', 'g2', 'h2']
trees_df['type'] = 'spe'
trees_df.loc[leaves, 'type'] = 'leaf'
#topo
topo_df = get_topo_time(trees_df)
full_topo_df = topo_df.join(trees_df)

# Wanted result:
assert (topo_df.topo_age[leaves] == 0).all()
assert (topo_df.topo_age[['c', 'b', 'd2', 'b2']] == 1).all()
assert (topo_df.topo_age[['a', 'a2']] == 2).all()
assert (topo_df.topo_age[['r', 'r2']] == 3).all()

assert (topo_df.topo_brlen[['e', 'f', 'g', 'h', 'c', 'a', 'e2', 'f2', 'g2', 'h2', 'd2', 'a2']] == 1).all()
assert (topo_df.topo_brlen[['d', 'b', 'c2', 'b2']] == 2).all()
assert (topo_df.topo_brlen[['r', 'r2']].isna().all())


topo_trees = parentdata_to_ete3(full_topo_df, 'topo_brlen')
