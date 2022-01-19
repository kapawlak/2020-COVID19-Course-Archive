#!/bin/bash
echo "Filtering sequences"
augur filter \
	--sequences data/sequences.fasta \
	--metadata data/metadata.tsv  \
	--output results/filtered.fasta \
	--min-date 2019 \
	&& \
echo "Aligning sequences"
augur align \
	--sequences results/filtered.fasta \
	--reference-sequence results/outgroup.fasta \
	--output results/aligned.fasta \
	&& \
echo "Building tree"
augur tree \
	--alignment results/aligned.fasta \
	--output results/tree_raw.nwk \
	&& \
augur refine --tree results/tree_raw.nwk \
	--alignment results/aligned.fasta \
	--metadata data/metadata.tsv \
	--output-tree results/tree.nwk \
	--output-node-data results/branch_lengths.json \
	--timetree \
	--coalescent opt \
	--date-inference marginal  \
	--clock-filter-iqd 4 \
	&& \
echo "Adding traits to tree"
augur traits \
	--tree results/tree.nwk \
	--metadata data/metadata.tsv \
	--output results/traits.json \
	--columns country \
	&& \
echo "Exporting results to auspice" 
augur export v2 \
	--tree results/tree.nwk \
	--metadata data/metadata.tsv \
	--node-data results/branch_lengths.json results/traits.json results/nt_muts.json \
	--output auspice/output.json \
	&& \
echo "View results on localhost:4000"
auspice view --datasetDir auspice &
