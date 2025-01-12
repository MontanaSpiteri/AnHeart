# Transcript niches

Code used for transcript-based niches adapted from Vannan, A., et al. Image-based spatial transcriptomics identifies molecular niche dysregulation associated with distal lung remodeling in pulmonary fibrosis. bioRxiv (2023).

1. Build mRNA spatial graphs per sample using radius threshold d. Each individual transcript is a node and edges are connections between these nodes (scripts: build_graph.py, build_graph_array.sh)
2. Create subgraph per sample from the randomly sample N root nodes (scripts: filter_and_subgraph.py, subgraph_array.sh)
3. Subgraphs from 9 AnHeart Xenium samples were merged into one graph (scripts: merge_subgraphs.py, merge_subgraphs_array.sh)
4. Two-hop GraphSAGE model was trained on the merged graph (scripts: train_subgraphs.py, train_subgraphs_job.sh)
5. Embedding prediction using trained model (scripts: embedd_graph.py, embedd_graph_job.sh)
6. pycave clustering to aggregate nodes for further transcript niche analysis (scripts: pycave_cluster.py)
7. Map niches back to transcripts (scripts: map_niches.R)
8. Aggregate transcript niches into hexbins (scripts: niche_hexbin.R)
9. Map nuclei to hexbins using a Euclidean distance (scripts: nuclei2hexbins.R)
