#!/bin/bash
mkdir -p data/ensemble_cdk2
for id in 6Q4G 6Q49 6Q4H 6Q48 6Q4J 6Q4E 6Q4K 6Q4D 6Q3B 6Q4I 6Q4B 6Q4A 6Q3F 6Q4F 4EK4 4FKL 2R3I 6Q3C 1GZ8; do
  code=$(curl -s -o "data/ensemble_cdk2/${id}.pdb" -w '%{http_code}' "https://files.rcsb.org/download/${id}.pdb")
  echo "${id}: http=${code}"
done
