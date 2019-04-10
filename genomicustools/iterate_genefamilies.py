#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from __future__ import print_function


from LibsDyogen import myTools, myProteinTree
import logging
logger = logging.getLogger(__name__)


if __name__ == '__main__':
    logging.basicConfig(format="%(levelname)s:%(message)s")
    args = myTools.checkArgs([('forestfile', myTools.File)], [], __doc__)

    for tree in myProteinTree.loadTree(args['forestfile']):
        rootinfo = tree.info[tree.root]
        try:
            tree_name = rootinfo['tree_name']
        except KeyError:
            logger.warning("No tree_name found at root: %d %s",
                           tree.root, rootinfo)
            tree_name = ';'.join((rootinfo.get('family_name', 'unnamed'),
                                  rootinfo['taxon_name']))
        print(tree_name)

