#!/bin/bash
for id in 1GZ8 4EK4 4FKL 6Q48 6Q49 6Q4A; do
  echo "=== ${id} ==="
  grep '^HETATM' "data/ensemble_cdk2/${id}.pdb" | awk '{print $4}' | sort -u
done
