"""Live progress display for `mvtpy build`.

Two modes, chosen automatically:

* **Interactive** (stdout is a terminal): a redrawing block - a progress bar,
  one line per unit currently running with its elapsed time, and a tally. The
  block is rewritten in place with ANSI cursor moves.
* **Plain** (piped, redirected, CI, ``TERM=dumb``, or ``MVT_NO_PROGRESS=1``):
  one durable line per event, exactly as before. Nothing redraws, nothing
  emits escape codes, and a log file stays greppable.

The plain mode is the important one to get right: pipelines here are usually
run under ``nohup``/``tee`` precisely because the machine has rebooted
mid-render before, and a progress widget that corrupts that log would be worse
than no progress display at all.

No third-party dependency - the port keeps to numpy/matplotlib, so this is
hand-rolled rather than pulling in rich or tqdm.
"""

from __future__ import annotations

import os
import shutil
import sys
import time
from typing import Dict, Optional, TextIO

__all__ = ["Progress", "format_duration"]

_SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
_FULL, _EMPTY = "█", "░"

#: Redraws per second in interactive mode. Fast enough to look live, slow
#: enough that a long build does not spend measurable time drawing.
_REDRAW_HZ = 8.0


def format_duration(seconds: float) -> str:
    """Compact, fixed-ish width: 45s, 6m12s, 1h04m."""
    seconds = max(0.0, float(seconds))
    if seconds < 60:
        return f"{seconds:.0f}s"
    if seconds < 3600:
        return f"{int(seconds // 60)}m{int(seconds % 60):02d}s"
    return f"{int(seconds // 3600)}h{int((seconds % 3600) // 60):02d}m"


def _interactive(stream: TextIO) -> bool:
    if os.environ.get("MVT_NO_PROGRESS"):
        return False
    if os.environ.get("TERM", "") in ("", "dumb"):
        return False
    try:
        return bool(stream.isatty())
    except (AttributeError, ValueError):
        return False


class Progress:
    """Progress reporter for a build run.

    The scheduler calls :meth:`started` / :meth:`finished` around each unit and
    :meth:`tick` while it waits; everything else is display.
    """

    def __init__(self, total: int, jobs: int = 1, stream: Optional[TextIO] = None,
                 interactive: Optional[bool] = None, fresh: int = 0):
        self.stream = stream or sys.stdout
        self.interactive = _interactive(self.stream) if interactive is None else interactive
        self.total = total
        self.jobs = jobs
        self.fresh = fresh
        self.started_at = time.time()
        self.running: Dict[str, tuple] = {}      # key -> (label, start)
        self.done = 0
        self.failed = 0
        self.durations: list = []
        self._lines_drawn = 0
        self._last_draw = 0.0

    # --- events -----------------------------------------------------------

    def info(self, message: str) -> None:
        """A durable line (kept above the live block in interactive mode)."""
        self._clear()
        self.stream.write(message + "\n")
        self.stream.flush()
        self._draw(force=True)

    def started(self, key: str, label: str, reason: str = "") -> None:
        self.running[key] = (label, time.time())
        if not self.interactive:
            self._plain(f"[mvt] build {label}" + (f": {reason}" if reason else ""))
        self._draw(force=True)

    def finished(self, key: str, label: str, ok: bool, elapsed: float,
                 detail: str = "") -> None:
        self.running.pop(key, None)
        self.done += 1
        if ok:
            self.durations.append(elapsed)
        else:
            self.failed += 1

        if self.interactive:
            self._clear()
            mark = "\x1b[32m✓\x1b[0m" if ok else "\x1b[31m✗\x1b[0m"
            self.stream.write(f"  {mark} {label}  {format_duration(elapsed)}\n")
            if detail:
                for line in detail.splitlines():
                    self.stream.write(f"      \x1b[2m{line}\x1b[0m\n")
            self.stream.flush()
            self._draw(force=True)
        else:
            if ok:
                self._plain(f"[mvt] done  {label} ({format_duration(elapsed)}) "
                            f"[{self.done}/{self.total}]")
            else:
                self._plain(f"[mvt] FAILED {label}")
                for line in detail.splitlines():
                    self._plain(f"       | {line}")

    def tick(self) -> None:
        """Called while the scheduler waits; refreshes elapsed times."""
        self._draw()

    def close(self, summary: str) -> None:
        self._clear()
        self.stream.write(summary + "\n")
        self.stream.flush()
        self._lines_drawn = 0

    # --- rendering --------------------------------------------------------

    def _plain(self, message: str) -> None:
        self.stream.write(message + "\n")
        self.stream.flush()

    def eta(self) -> Optional[float]:
        """Seconds remaining, from mean unit duration spread over the workers."""
        if not self.durations or self.done >= self.total:
            return None
        mean = sum(self.durations) / len(self.durations)
        remaining = self.total - self.done
        return mean * remaining / max(1, self.jobs)

    def bar(self, width: int = 30) -> str:
        fraction = self.done / self.total if self.total else 1.0
        filled = int(round(fraction * width))
        return _FULL * filled + _EMPTY * (width - filled)

    def _render(self) -> list:
        percent = (self.done / self.total * 100) if self.total else 100.0
        elapsed = time.time() - self.started_at
        eta = self.eta()
        header = (f"[{self.bar()}] {self.done}/{self.total} {percent:3.0f}%  "
                  f"{format_duration(elapsed)} elapsed")
        if eta is not None:
            header += f"  ~{format_duration(eta)} left"

        frame = _SPINNER[int(time.time() * 10) % len(_SPINNER)]
        lines = [header]
        for label, start in sorted(self.running.values(), key=lambda item: item[1]):
            lines.append(f"  \x1b[36m{frame}\x1b[0m {label}  "
                         f"\x1b[2m{format_duration(time.time() - start)}\x1b[0m")

        tally = f"  {self.done - self.failed} built"
        if self.failed:
            tally += f", \x1b[31m{self.failed} failed\x1b[0m"
        if self.fresh:
            tally += f", {self.fresh} up to date"
        queued = self.total - self.done - len(self.running)
        if queued > 0:
            tally += f", {queued} queued"
        lines.append(tally)

        width = shutil.get_terminal_size((100, 24)).columns
        return [line[:width + 20] for line in lines]   # +20 for escape sequences

    def _draw(self, force: bool = False) -> None:
        if not self.interactive:
            return
        now = time.time()
        if not force and now - self._last_draw < 1.0 / _REDRAW_HZ:
            return
        self._last_draw = now
        self._clear()
        lines = self._render()
        self.stream.write("\n".join(lines) + "\n")
        self.stream.flush()
        self._lines_drawn = len(lines)

    def _clear(self) -> None:
        """Erase the live block so the next write lands on a clean area."""
        if not self.interactive or not self._lines_drawn:
            return
        self.stream.write(f"\x1b[{self._lines_drawn}A")   # up N lines
        self.stream.write("\x1b[J")                        # clear to end of screen
        self.stream.flush()
        self._lines_drawn = 0
