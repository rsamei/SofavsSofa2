"""Run the whole pipeline in order: build the upstream tables, the hourly scores and the cohort, then the analyses.

Usage: BQ_PROJECT=<project> python run_all.py [--from N] [--skip-existing]
  --from N         start at step N, e.g. to resume after an interruption;
  --skip-existing  in step 1, keep upstream tables that already exist.
Results are written as JSON to OUTPUT_DIR (default output/).
"""
import subprocess
import sys
from pathlib import Path

STEPS = sorted((Path(__file__).resolve().parent / "pipeline").glob("[0-9][0-9]_*.py"))


def main():
    start = int(sys.argv[sys.argv.index("--from") + 1]) if "--from" in sys.argv else 1
    for script in STEPS:
        n = int(script.name[:2])
        if n < start:
            continue
        args = [sys.executable, str(script)] + (["--skip-existing"] if n == 1 and "--skip-existing" in sys.argv else [])
        print(f"== step {n}: {script.name}", flush=True)
        if subprocess.run(args).returncode != 0:
            print(f"step {n} failed; rerun with --from {n}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
