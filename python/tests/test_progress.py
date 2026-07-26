"""Tests for the build progress display.

The behaviour that matters most is the plain (non-terminal) mode: builds here
are routinely run under nohup/tee because the machine has rebooted mid-run
before, and a progress widget that sprayed escape codes into that log would be
worse than no display at all.
"""

from __future__ import annotations

import io
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy.progress import Progress, format_duration  # noqa: E402


class FakeTTY(io.StringIO):
    def isatty(self):
        return True


def test_format_duration_units():
    assert format_duration(0) == "0s"
    assert format_duration(45) == "45s"
    assert format_duration(372) == "6m12s"
    assert format_duration(3900) == "1h05m"
    assert format_duration(-5) == "0s"          # clock skew must not print "-5s"


# --- plain mode -----------------------------------------------------------


def test_plain_mode_emits_no_escape_codes():
    stream = io.StringIO()
    progress = Progress(total=2, jobs=2, stream=stream, interactive=False)
    progress.started("a", "slim #00", reason="missing output")
    progress.finished("a", "slim #00", ok=True, elapsed=97.0)
    progress.started("b", "slim #01")
    progress.finished("b", "slim #01", ok=False, elapsed=3.0, detail="Traceback\nBoom")
    progress.close("[mvt] 1 built, 1 failed")

    output = stream.getvalue()
    assert "\x1b" not in output, "escape codes leaked into a non-terminal stream"
    assert "[mvt] build slim #00: missing output" in output
    assert "[mvt] done  slim #00 (1m37s) [1/2]" in output
    assert "[mvt] FAILED slim #01" in output
    assert "       | Boom" in output


def test_plain_mode_never_redraws():
    """Every line must be durable: no cursor moves, one line per event."""
    stream = io.StringIO()
    progress = Progress(total=1, jobs=1, stream=stream, interactive=False)
    progress.started("a", "gps")
    for _ in range(5):
        progress.tick()
    progress.finished("a", "gps", ok=True, elapsed=1.0)
    assert stream.getvalue().count("\n") == 2       # exactly the two events


def test_environment_can_force_plain_mode(monkeypatch):
    monkeypatch.setenv("MVT_NO_PROGRESS", "1")
    monkeypatch.setenv("TERM", "xterm-256color")
    assert not Progress(total=1, stream=FakeTTY()).interactive


def test_dumb_terminal_is_plain(monkeypatch):
    monkeypatch.delenv("MVT_NO_PROGRESS", raising=False)
    monkeypatch.setenv("TERM", "dumb")
    assert not Progress(total=1, stream=FakeTTY()).interactive


def test_non_tty_is_plain(monkeypatch):
    monkeypatch.delenv("MVT_NO_PROGRESS", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    assert not Progress(total=1, stream=io.StringIO()).interactive


# --- interactive mode -----------------------------------------------------


def interactive_progress(total=4, jobs=2):
    return Progress(total=total, jobs=jobs, stream=FakeTTY(), interactive=True)


def test_interactive_draws_a_bar_and_running_units():
    progress = interactive_progress()
    progress.started("a", "slim 2022-11-18 #00")
    progress.started("b", "slim 2022-11-18 #01")
    output = progress.stream.getvalue()

    assert "█" in output or "░" in output
    assert "slim 2022-11-18 #00" in output
    assert "0/4" in output


def test_bar_fills_with_completion():
    progress = interactive_progress(total=4)
    assert progress.bar(width=4) == "░░░░"
    progress.done = 2
    assert progress.bar(width=4) == "██░░"
    progress.done = 4
    assert progress.bar(width=4) == "████"


def test_completed_units_scroll_above_the_live_block():
    progress = interactive_progress()
    progress.started("a", "gps 2022-11-18")
    progress.finished("a", "gps 2022-11-18", ok=True, elapsed=700.0)
    output = progress.stream.getvalue()
    assert "✓" in output and "gps 2022-11-18" in output and "11m40s" in output


def test_failed_unit_shows_detail():
    progress = interactive_progress()
    progress.started("a", "slim #03")
    progress.finished("a", "slim #03", ok=False, elapsed=2.0, detail="ValueError: bad")
    assert "✗" in progress.stream.getvalue()
    assert "ValueError: bad" in progress.stream.getvalue()


def test_eta_uses_observed_durations_over_workers():
    progress = interactive_progress(total=10, jobs=2)
    assert progress.eta() is None                  # nothing measured yet
    for index in range(2):
        progress.finished(f"u{index}", f"unit {index}", ok=True, elapsed=100.0)
    # 8 remaining x 100s, spread over 2 workers.
    assert progress.eta() == pytest.approx(400.0)


def test_eta_is_none_once_complete():
    progress = interactive_progress(total=1, jobs=1)
    progress.finished("a", "unit", ok=True, elapsed=10.0)
    assert progress.eta() is None


def test_info_lines_survive_the_redraw():
    progress = interactive_progress()
    progress.started("a", "slim #00")
    progress.info("[mvt] up to date: fields 2022-11-18")
    assert "up to date: fields" in progress.stream.getvalue()


def test_close_clears_the_block_and_prints_the_summary():
    progress = interactive_progress()
    progress.started("a", "slim #00")
    progress.close("[mvt] 1 built, 0 failed")
    assert progress.stream.getvalue().endswith("[mvt] 1 built, 0 failed\n")
