"""Ensemble analysis of docked pose collections."""

from __future__ import annotations

from typing import Iterable, List


def cluster_poses(poses: Iterable, rmsd_cutoff: float = 2.0) -> List[List]:
    """Cluster a set of poses by structural similarity (e.g. RMSD).

    Args:
        poses: An iterable of pose objects.
        rmsd_cutoff: RMSD threshold (in Angstroms) for grouping poses.

    Returns:
        A list of clusters, each a list of poses.
    """
    raise NotImplementedError


def rank_poses(poses: Iterable, scoring_key: str = "score") -> List:
    """Rank poses by a scoring metric.

    Args:
        poses: An iterable of pose objects.
        scoring_key: The attribute/key to rank poses by.

    Returns:
        Poses sorted best-to-worst by the scoring key.
    """
    raise NotImplementedError


def compute_consensus_pose(poses: Iterable):
    """Compute a consensus/representative pose from an ensemble.

    Args:
        poses: An iterable of pose objects.

    Returns:
        A single representative pose.
    """
    raise NotImplementedError
