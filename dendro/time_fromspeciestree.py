#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Set up branch lengths using species tree branch lengths"""


from sys import stdout
import argparse as ap
import ete3
from LibsDyogen import myPhylTree
from dendro.parsers import read_multinewick
#from seqtools.specify import tree_specify  
import logging
logger = logging.getLogger(__name__)

ENSEMBL_VERSION = 93


#~~> See prune2family and generate_dNdS
def time_fromspeciestree(tree, phyltree):
    """Assume the tree has a 'S' feature storing the species name"""
    
    for node in tree.traverse('preorder'):
        #event = get_event(node)
        #print(node.name, node.features)
        if node.D == 'N' and not node.features.intersection(('A', 'age')):
            node_species = node.S.replace('.', ' ')
            node.add_feature('age', phyltree.ages[node_species])

            logger.info('* %s: S=%s; D=%s; age=%s', node.name, node.S, node.D, node.age)
            if not node.is_root():
                parent = node.up

                intermediates = [[node]]  # One list per taxon lineage
                inter_species = [node_species]
                while parent is not None and parent.D == 'Y':
                    parent_species = parent.S.replace('.', ' ')
                    if parent_species == node_species:
                        intermediates[-1].append(parent)
                    else:
                        logger.warning("Intermediate speciations are missing: "
                                       "%s -> %s",
                                       parent.name, node.name)
                        #node = parent
                        inter_species.append(parent_species)
                        intermediates.append([parent])  # New taxon lineage
                        node_species = parent_species

                    parent = parent.up

                logger.debug('%d %s: %d %s', len(inter_species), inter_species,
                             len(intermediates), [[n.name for n in interl]
                                                  for interl in intermediates])
                inter_species.append((phyltree.parent[node_species].name
                                      if parent is None else
                                      parent.S.replace('.', ' ')))

                node_species = inter_species.pop(0)
                for k, (parent_species, inter_lineage) in enumerate(zip(inter_species,
                                                                  intermediates)):
                    # process intermediates
                    step = phyltree.ages[parent_species] - phyltree.ages[node_species]
                    if k == 0:
                        step /= len(inter_lineage)
                    else:
                        # The starting speciation node is not included.
                        # There are only dup nodes between ghost speciations:
                        # Divide the interval length by N+1.
                        step /= (len(inter_lineage) + 1)
                        # Also add this step to the node in the previous lineage.
                        intermediates[k-1][-1].dist += step
                    logger.debug('  > k=%d %s->%s step=(%s-%s)/%s=%s', k,
                                 parent_species, node_species,
                                 phyltree.ages[parent_species],
                                 phyltree.ages[node_species],
                                 len(inter_lineage) if k==0 else (len(inter_lineage)+1),
                                 step)

                    for n in inter_lineage:
                        n.dist = step

                    node_species = parent_species


def time_fromspeciestreeIO(treefile, phyltreefile, outfile=None):
    phyltree = myPhylTree.PhylogeneticTree(phyltreefile)
    with open(treefile) as f:
        lines = f.readlines()
    with (open(outfile, 'w') if outfile else stdout) as out:
        for treetxt in read_multinewick(lines):
            tree = ete3.Tree(treetxt, format=1)
            time_fromspeciestree(tree, phyltree)
            newick = tree.write(format=1, format_root_node=False)
            out.write(newick + '\n')


if __name__ == '__main__':
    logging.basicConfig(format=logging.BASIC_FORMAT)
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('phyltreefile')
    parser.add_argument('treefile')
    parser.add_argument('outfile', nargs='?')
    parser.add_argument('-v', '--verbose', action='count', default=0)
    
    dargs = vars(parser.parse_args())
    verbosity = dargs.pop('verbose')
    if verbosity > 1:
        logger.setLevel(logging.DEBUG)
    elif verbosity > 0:
        logger.setLevel(logging.INFO)

    time_fromspeciestreeIO(**dargs)
