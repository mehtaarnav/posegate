# posegate/tests/test_run_obabel.py
import subprocess
import sys
import os

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from prep_ensemble import run_obabel


def test_run_obabel_gives_actionable_error_when_command_not_found(monkeypatch):
    """Found via an actual clean-environment install test (pip install .
    alone does not provide obabel, since it has no bundled pip wheel):
    the bare subprocess error ('exit status 127') gave no indication
    obabel was even the problem. Must raise a message naming obabel and
    how to install it, not just propagate the raw exit code."""
    def fake_run(*a, **k):
        raise subprocess.CalledProcessError(
            returncode=127, cmd='obabel ...', stderr=b'')
    monkeypatch.setattr(subprocess, 'run', fake_run)

    with pytest.raises(RuntimeError, match="obabel"):
        run_obabel('in.pdb', 'out.sdf')


def test_run_obabel_other_failure_includes_stderr(monkeypatch):
    def fake_run(*a, **k):
        raise subprocess.CalledProcessError(
            returncode=1, cmd='obabel ...', stderr=b'some real obabel error')
    monkeypatch.setattr(subprocess, 'run', fake_run)

    with pytest.raises(RuntimeError, match="some real obabel error"):
        run_obabel('in.pdb', 'out.sdf')
