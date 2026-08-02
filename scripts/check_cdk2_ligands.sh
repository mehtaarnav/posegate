#!/bin/bash
for f in data/ensemble_cdk2/*.pdb; do
  id=$(basename "$f" .pdb)
  lig=$(grep '^HETATM' "$f" | awk '{print $4}' | sort -u | grep -Ev '^(HOH|SO4|GOL|EDO|DMS|PEG|PO4|CL|NA|K|MG|CA|ZN|ACT|TRS|BME|IMD|MPD|1PE|PGE|MRD|FMT|EPE)$' | tr '\n' ',')
  echo "${id}: ${lig}"
done
