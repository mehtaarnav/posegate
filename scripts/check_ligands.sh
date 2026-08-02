#!/bin/bash
for f in data/ensemble/*.pdb data/3mxf.pdb; do
  echo "=== ${f} ==="
  grep '^HETATM' "${f}" | awk '{print $4}' | sort -u | grep -Ev '^(HOH|SO4|GOL|EDO|DMS|PEG|PO4|CL|NA|K|MG|CA|ZN|ACT|TRS|BME|IMD|MPD|1PE|PGE)$'
done
