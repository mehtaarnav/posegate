"""Shared geometry and I/O helper utilities."""

from __future__ import annotations

from typing import Sequence


def load_structure(path: str):
    """Load a molecular structure from disk (PDB, SDF, MOL2, etc.).

    Args:
        path: Path to the structure file.

    Returns:
        A parsed molecule/structure object.
    """
    raise NotImplementedError


def euclidean_distance(point_a: Sequence[float], point_b: Sequence[float]) -> float:
    """Compute the Euclidean distance between two 3D points.

    Args:
        point_a: (x, y, z) coordinates of the first point.
        point_b: (x, y, z) coordinates of the second point.

    Returns:
        The distance between the two points.
    """
    raise NotImplementedError


def compute_angle(point_a: Sequence[float], point_b: Sequence[float], point_c: Sequence[float]) -> float:
    """Compute the angle (in degrees) at vertex ``point_b`` formed by A-B-C.

    Args:
        point_a: Coordinates of the first point.
        point_b: Coordinates of the vertex point.
        point_c: Coordinates of the third point.

    Returns:
        The angle in degrees.
    """
    raise NotImplementedError
