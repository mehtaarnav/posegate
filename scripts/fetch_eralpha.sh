#!/bin/bash
for id in 6VPF 2AYR 7KBS; do
  curl -s -o "data/five_family/${id}.pdb" "https://files.rcsb.org/download/${id}.pdb"
  echo "=== ${id} ==="
  grep '^HETATM' "data/five_family/${id}.pdb" | awk '{print $4}' | sort -u | grep -Ev '^(HOH|SO4|GOL|EDO|DMS|PEG|PO4|CL|NA|K|MG|CA|ZN|ACT|TRS|BME|IMD|MPD|1PE|PGE|MRD|FMT|EPE)$'
done
