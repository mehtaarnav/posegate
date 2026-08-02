#!/bin/bash
grep '^HETATM' data/ensemble_cdk2/3EZR.pdb | grep 'EZR' | awk '{print $5, $6}' | sort -u
