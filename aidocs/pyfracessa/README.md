# PyFracESSA Timing Reference

Status: current internal benchmark workflow.

Public installation, matrix input, result fields, multiprocessing, and sink behavior are documented in `../../README.md` and the
generated pages under `../../docs/`. This file intentionally does not duplicate that user documentation.

## Purpose

`python/pyfracessa/timing.py` measures one build per invocation, one matrix at a time, and stores native analyzer durations in the
canonical SQLite database. It does not call `run_multiprocessing()`. Reuse one session name across invocations to compare builds
without loading incompatible native extensions into the same Python process.

Canonical local comparisons use:

- CPU 2;
- the persistent Pybind backend;
- a freshly built binary for each revision;
- stored per-matrix calibrations;
- a 0.5-second target per matrix and method unless the experiment records another target; and
- the median native `elapsed_ns`, not Python wall time.

Never compare persistent-Pybind timing with CLI process timing.

## Canonical run

From the repository root:

```bash
PYTHONPATH=python python3 -m pyfracessa.timing run \
  --backend pybind --module-dir cpp/build \
  --build-label main --source-ref main --revision "$(git rev-parse HEAD)" \
  --method fast --method safe \
  --cpu 2 --comment "before candidate change"
```

The command prints its generated session. Pass that same session with a different build label when measuring the comparison build.
Then report the session:

```bash
PYTHONPATH=python python3 -m pyfracessa.timing report latest --baseline main
```

Use the actual printed session instead of `latest` when concurrent or later timing runs could make `latest` ambiguous.

## Matrix selection and calibration

The default selection is the `small` size class. `--size-class` also accepts `medium`, `large`, `super_large`, and `all`; repeated
`--matrix-id` arguments select exact database rows.

Each matrix and method starts from `fast_calibration_ns` or `safe_calibration_ns` in `testdata/fracessa_testdata.sqlite3`. A positive
calibration chooses `ceil(target / calibration)` samples. A calibration at or above the target, or the `-1` timeout sentinel, chooses
one sample. Missing calibration data is an error. Change the default 0.5-second target with `--target-seconds` and record that
deviation in the benchmark report.

The native module remains loaded for every selected method and matrix in one invocation. Python wall time is stored as metadata but
does not select the sample count or reported result.

## Stored provenance and checks

Every timing row records the moving `source_ref`, exact `revision`, module or executable SHA-256, backend, build label, machine, CPU,
target, iteration count, native median, measured wall duration, and safe-fallback reason. The runner compares the observed ESS count
with the canonical database count and leaves unsafe mismatches visible rather than discarding them.

Reports include matrix dimension, circular-symmetry status, and the lower bound
`gamma_lower_bound = expected_ess ** (1 / dimension)`. Ratio reports compare only matrices shared with the named baseline.

## Legacy CLI measurements

Current CLI output reports nanoseconds directly. Historical executables may need an explicit unit and default-method declaration:

- use `--fast-default --cli-unit s` when the old no-flag route corresponds to current `fast` and reports seconds;
- use `--safe-default` when the old no-flag route corresponds to current `safe`.

CLI runs start a child process for every sample. Keep them in separate reports from persistent-Pybind microbenchmarks. The adapter
maps current `fast` and `safe` labels onto the supported historical mode or Boolean interfaces.

## Command reference

Use the implementation's help output for the complete current option list:

```bash
PYTHONPATH=python python3 -m pyfracessa.timing --help
PYTHONPATH=python python3 -m pyfracessa.timing run --help
PYTHONPATH=python python3 -m pyfracessa.timing report --help
```
