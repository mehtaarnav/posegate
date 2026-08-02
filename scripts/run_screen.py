#!/usr/bin/env python
"""Entry-point script for running a docking screen.

Usage:
    python scripts/run_screen.py --receptor receptor.pdb --ligands ligands/ --out results/
"""

from __future__ import annotations

import argparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a docking screen with posegate.")
    parser.add_argument("--receptor", required=True, help="Path to the receptor structure.")
    parser.add_argument("--ligands", required=True, help="Path to a ligand file or directory.")
    parser.add_argument("--out", required=True, help="Output directory for results.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raise NotImplementedError


if __name__ == "__main__":
    main()
