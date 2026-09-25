"""Shared stdout logger for the LocusBlend package.

Single definition of the log line format used by app.py and the
locusblend_web modules.
"""

import time


def log(msg):
    """Print a timestamped message to stdout."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)
