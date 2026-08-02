#!/bin/bash
# CDK2 ensemble for the conserved-contact miner and the visGReMLIN
# comparison (scripts/compare_visgremlin.py). 22 structures.
#
# The composition of this list matters more than its length. 6Q3B-6Q4K
# are a single fragment-screen deposition series that binds the
# hinge/adenine subpocket without reaching the DFG region. Mining the
# 19-structure list this script previously fetched, which was dominated
# by that series, dropped ASP145 from the miner's output entirely, and
# restricting the series to its drug-like members (>=10 heavy atoms) did
# not recover it. The last three IDs are inhibitors that do reach the
# residues the 6Q series misses, and they restore reference-residue
# coverage to 9/9, so removing them to make the ensemble more uniform
# would undo that.
mkdir -p data/ensemble_cdk2
for id in 6Q4G 6Q49 6Q4H 6Q48 6Q4J 6Q4E 6Q4K 6Q4D 6Q3B 6Q4I 6Q4B 6Q4A 6Q3F 6Q4F 4EK4 4FKL 2R3I 6Q3C 1GZ8 5JQ5 3EZR 3WBL; do
  code=$(curl -s -o "data/ensemble_cdk2/${id}.pdb" -w '%{http_code}' "https://files.rcsb.org/download/${id}.pdb")
  echo "${id}: http=${code}"
done
