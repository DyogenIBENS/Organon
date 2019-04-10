/*
Query the list of (node ids, tree ids) corresponding to a gene_split
===================================================================

Run with:
mysql -h ensembldb.ensembl.org -u anonymous ensembl_compara_85 <script.sql >output.tab

# Description of the useful Ensembl database tables and entries ##

See http://jul2016.archive.ensembl.org/info/docs/api/compara/compara_schema.html

homology table:
 - description "gene_split"
 - gene_tree_node_id
 - gene_tree_root_id

Or gene_tree_node_attr (gene split only once)
 - node_id
 - node_type "gene_split"

gene_tree_node table
 - node_id
 - root_id
 - left_index/right_index (what about multifurcations)
 - seq_member_id

gene_tree_root table
 - root_id
 - stable_id
 - (gene_align_id)
*/

SELECT n.node_id, n.root_id, stable_id, node_type, member_type, tree_type, n.left_index, n.right_index, s.node_name AS species_name
FROM gene_tree_node            AS n
INNER JOIN gene_tree_root      AS r ON (n.root_id = r.root_id)
INNER JOIN gene_tree_node_attr AS a ON (n.node_id = a.node_id)
LEFT JOIN species_tree_node    AS s ON (a.species_tree_node_id = s.node_id)
WHERE (node_type="gene_split"
   AND member_type="protein"
   AND tree_type="tree"
   AND stable_id IS NOT NULL)
ORDER BY n.root_id;
