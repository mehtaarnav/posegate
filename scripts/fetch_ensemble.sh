#!/bin/bash
set -e
mkdir -p data/ensemble
for id in 6I7Y 8YMG 8YMI 7ZE6 4ZC9; do
  code=$(curl -s -o "data/ensemble/${id}.pdb" -w '%{http_code}' "https://files.rcsb.org/download/${id}.pdb")
  n=$(wc -l < "data/ensemble/${id}.pdb")
  echo "${id}: http=${code} lines=${n}"
done
