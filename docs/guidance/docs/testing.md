# Regression and Smoke Tests

The `tests/` directory contains ready-to-run cases and helper scripts for quick
health checks and numerical regression checks.

Use smoke tests for rapid validation after build or input-generation changes. Use regression tests before merging solver, discretization, I/O, or physics modifications.

## Test Scripts

| Script | Purpose |
| --- | --- |
| `tests/run_smoke.sh` | Execute a brief 10-iteration subset of representative cases for rapid validation |
| `tests/run_regression.sh` | Run or validate standard/extended regression suites against reference metrics |
| `tests/update_reference.sh` | Replace `reference.json` files with newly generated metrics (use only after validated changes) |
| `tests/tools/check_metrics.py` | Compare `regression_test_metrics.json` against `reference.json` using specified tolerances |
| `tests/tools/tolerances.json` | Metric comparison tolerance specifications |

## Smoke Tests

From the `tests/` directory:

```bash
cd tests
bash run_smoke.sh
```

The smoke suite is designed to detect build failures, missing files, invalid input generation, and early runtime errors. It is not intended as a substitute for statistically converged production validation.

## Regression Tests

From the `tests/` directory, execute:

```bash
cd tests
bash run_regression.sh
```

The script supports two operational modes:

| Mode | Behavior |
| --- | --- |
| `run` | Execute each case, await `regression_test_metrics.json`, then validate metrics |
| `check` | Skip solver execution and validate existing metrics only |

Two test suites are available:

| Suite | Scope |
| --- | --- |
| `standard` | Compact default set for routine validation |
| `extended` | Comprehensive set covering thermal, periodic, inlet/outlet, annular, and pipe cases |

In continuous integration environments, defaults are applied without interactive user prompts.

## Reference Metrics

Each regression case should contain:

| File | Meaning |
| --- | --- |
| `input_chapsim.ini` | Case setup. |
| `run_chapsim.sh` | Case run wrapper. |
| `reference.json` | Accepted reference metrics. |
| `regression_test_metrics.json` | Metrics produced by the latest run. |

The comparison is performed by:

```bash
python3 tests/tools/check_metrics.py \
  tests/<case>/regression_test_metrics.json \
  tests/<case>/reference.json \
  tests/tools/tolerances.json
```

## Updating References

Update references only after confirming that any numerical changes are intentional and physically validated.

```bash
cd tests
bash update_reference.sh
```

This operation copies each available `regression_test_metrics.json` to `reference.json`. Review the resulting version-control diff before committing.

## Recommended Testing Policy

- **Input-generation tool modifications**: Execute smoke tests and at least the affected generated cases
- **I/O, restart, mesh, or boundary-condition modifications**: Execute smoke tests plus the standard regression suite
- **Numerical scheme or physics model modifications**: Execute the extended regression suite and inspect representative monitor histories
- **Reference metric updates**: Document the rationale explaining why new metrics are expected
